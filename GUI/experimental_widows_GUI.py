import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import subprocess
import os
import threading
import shlex
import re

# List to store directories for mounting
mount_directories = []
docker_process = None  # Global variable to track the running Docker process

# Lists to store selected BAM files
test_samples_files = []
baseline_samples_files = []

def windows_to_unix_path(path):
    """Convert a Windows path to a Unix-style path for Docker."""
    path = os.path.abspath(path)
    path = path.replace('\\', '/')  # Convert backslashes to forward slashes
    
    # Handle drive letters (R: -> /mnt/r)
    match = re.match(r'^([a-zA-Z]):', path)
    if match:
        drive_letter = match.group(1).lower()  # Extract drive letter and convert to lowercase
        unix_path = f'/mnt/{drive_letter}{path[2:]}'  # Replace 'C:' with '/mnt/c' and similar
    elif path.startswith('//') or path.startswith('\\\\'):
        # Handle UNC paths (e.g., \\server\share)
        path = path.lstrip('/').lstrip('\\')
        unix_path = f'/mnt/network/{path}'
    else:
        # Handle relative paths
        unix_path = path  # Relative paths do not need conversion
    
    return unix_path

def normalize_path(path):
    """Normalize user-supplied paths to handle mixed slashes."""
    return path.replace('\\', '/')

def browse_file(entry):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_path = normalize_path(file_path)
        entry.delete(0, tk.END)
        entry.insert(0, file_path)
        add_to_mount_list(file_path)

def browse_directory(entry):
    folder_path = filedialog.askdirectory()
    if folder_path:
        folder_path = normalize_path(folder_path)
        entry.delete(0, tk.END)
        entry.insert(0, folder_path)
        add_to_mount_list(folder_path)

def select_test_samples():
    global test_samples_files
    files = filedialog.askopenfilenames(filetypes=[("BAM files", "*.bam")])
    if files:
        test_samples_files = list(files)
        # Update the label to show the number of files selected
        test_samples_label.config(text=f"{len(test_samples_files)} files selected")
        # Add directories to mount list
        for file_path in test_samples_files:
            add_to_mount_list(file_path)

def select_baseline_samples():
    global baseline_samples_files
    files = filedialog.askopenfilenames(filetypes=[("BAM files", "*.bam")])
    if files:
        baseline_samples_files = list(files)
        # Update the label to show the number of files selected
        baseline_samples_label.config(text=f"{len(baseline_samples_files)} files selected")
        # Add directories to mount list
        for file_path in baseline_samples_files:
            add_to_mount_list(file_path)

def add_to_mount_list(file_path):
    """Extract the directory of the file and add it to the list if not already present."""
    file_path = normalize_path(file_path)
    directory = os.path.dirname(file_path)
    if directory not in mount_directories:
        mount_directories.append(directory)

def build_docker_volumes():
    """Build the Docker volume list with proper quoting and unique container paths."""
    volume_mounts = []
    for directory in mount_directories:
        unix_path = windows_to_unix_path(directory)
        volume_mounts.extend(["-v", f"{directory}:{unix_path}"])
    return volume_mounts

def stream_output(process):
    """Stream output from the subprocess and display it in real-time."""
    console_output.configure(state='normal')
    for line in iter(process.stdout.readline, b''):
        console_output.insert(tk.END, line.decode('utf-8'))
        console_output.see(tk.END)  # Auto-scroll to the end
    process.stdout.close()
    console_output.configure(state='disabled')
    
    # Update status after the process completes
    container_status_label.config(text="Finished", foreground="green")
    status_label.config(text="Completed")

def validate_inputs():
    """Validate that all required inputs are provided and exist."""
    global targets, annotation, output_dir

    # Normalize paths
    output_dir = normalize_path(output_dir_entry.get())
    targets = normalize_path(targets_entry.get()) if targets_entry.get() else "/data/default_targets.bed"
    annotation = normalize_path(annotation_entry.get()) if annotation_entry.get() else "/data/default_annotation.gff"

    if not output_dir:
        messagebox.showerror("Error", "Output Directory is required.")
        return False
    if not os.path.exists(output_dir):
        messagebox.showerror("Error", f"Output Directory does not exist:\n{output_dir}")
        return False

    # Check optional files
    if targets_entry.get() and not os.path.exists(targets):
        messagebox.showerror("Error", f"Target File does not exist:\n{targets}")
        return False
    if annotation_entry.get() and not os.path.exists(annotation):
        messagebox.showerror("Error", f"Annotation File does not exist:\n{annotation}")
        return False

    # Check that test_samples_files and baseline_samples_files are not empty
    if not test_samples_files:
        messagebox.showerror("Error", "No test samples selected.")
        return False

    if not baseline_samples_files:
        messagebox.showerror("Error", "No baseline samples selected.")
        return False

    return True

def run_docker():
    global docker_process  # Access and modify the global docker_process variable

    if not validate_inputs():
        return

    if not docker_installed:
        messagebox.showerror("Error", "Docker is not installed.")
        return

    if not container_downloaded:
        messagebox.showerror("Error", "Docker container is not downloaded.")
        return

    # At this point, targets, annotation, output_dir are set and normalized
    # test_samples_files and baseline_samples_files are set

    # Add directories to mount list
    add_to_mount_list(output_dir)
    if targets_entry.get():
        add_to_mount_list(targets_entry.get())
    if annotation_entry.get():
        add_to_mount_list(annotation_entry.get())

    # Create test samples file and baseline samples file in the output directory
    test_samples_file_path = os.path.join(output_dir, 'test_samples.txt')
    baseline_samples_file_path = os.path.join(output_dir, 'baseline_samples.txt')

    # Convert .bam file paths to Docker paths and write to files
    try:
        with open(test_samples_file_path, 'w') as test_file:
            for bam_path in test_samples_files:
                docker_path = windows_to_unix_path(bam_path)
                test_file.write(docker_path + '\n')

        with open(baseline_samples_file_path, 'w') as baseline_file:
            for bam_path in baseline_samples_files:
                docker_path = windows_to_unix_path(bam_path)
                baseline_file.write(docker_path + '\n')
    except Exception as e:
        messagebox.showerror("Error", f"Failed to create sample files: {e}")
        return

    # Build the docker volumes list
    docker_volumes = build_docker_volumes()

    # Construct the Docker command as a list
    docker_command = [
        "docker", "run", *docker_volumes,
        "-it", "murphydaviducl/runexomedepth:latest",
        "Rscript", "/runExomeDepth/ExomeDepth.R",
        "--targets", windows_to_unix_path(targets),
        "--annotation", windows_to_unix_path(annotation),
        "--test-samples", windows_to_unix_path(test_samples_file_path),
        "--baseline-samples", windows_to_unix_path(baseline_samples_file_path),
        "--output-directory", windows_to_unix_path(output_dir)
    ]

    console_output.configure(state='normal')
    console_output.insert(tk.END, f"Running Docker Command:\n{' '.join(shlex.quote(arg) for arg in docker_command)}\n\n")
    console_output.see(tk.END)
    console_output.configure(state='disabled')

    status_label.config(text="Running Docker...")
    container_status_label.config(text="Running", foreground="orange")

    try:
        # Start the subprocess with stdout and stderr piped
        docker_process = subprocess.Popen(docker_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
        
        # Create a thread to stream the output from the process to the console_output widget
        thread = threading.Thread(target=stream_output, args=(docker_process,))
        thread.start()

    except Exception as e:
        container_status_label.config(text="Error", foreground="red")
        console_output.configure(state='normal')
        console_output.insert(tk.END, f"Error: {str(e)}\n")
        console_output.see(tk.END)
        console_output.configure(state='disabled')

def on_closing():
    """Handle the window closing event to terminate the Docker process if running."""
    global docker_process

    if docker_process is not None and docker_process.poll() is None:  # Check if the process is running
        if messagebox.askokcancel("Quit", "The Docker process is still running. Do you want to terminate it and exit?"):
            docker_process.terminate()  # Send SIGTERM
            try:
                docker_process.wait(timeout=5)  # Give the process some time to terminate
            except subprocess.TimeoutExpired:
                docker_process.kill()  # Forcefully kill if it doesn't terminate in time
        else:
            return  # Don't close the window if the user cancels the quit action

    root.destroy()  # Close the application

def check_docker_installed():
    """Check if Docker is installed by running 'docker --version'."""
    try:
        subprocess.run(["docker", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        docker_status_label.config(text="Docker Installed", foreground="green")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        docker_status_label.config(text="Docker Not Installed", foreground="red")
        return False

def check_container_downloaded():
    """Check if the required Docker container is downloaded."""
    try:
        result = subprocess.run(["docker", "images", "-q", "murphydaviducl/runexomedepth:latest"], 
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.stdout.strip():
            container_status_label.config(text="Ready to Run", foreground="green")
            return True
        else:
            container_status_label.config(text="Not Downloaded", foreground="red")
            return False
    except subprocess.CalledProcessError:
        container_status_label.config(text="Error Checking Container", foreground="red")
        return False

def download_container():
    """Download the Docker container if it's not downloaded."""
    if not docker_installed:
        messagebox.showerror("Error", "Docker is not installed.")
        return

    container_status_label.config(text="Downloading...", foreground="orange")
    status_label.config(text="Downloading Container...")
    
    def stream_output(process):
        """Stream the output from the process to the console_output box."""
        console_output.configure(state='normal')
        for line in iter(process.stdout.readline, b''):
            console_output.insert(tk.END, line.decode('utf-8'))
            console_output.see(tk.END)  # Auto-scroll to the end
        process.stdout.close()
        console_output.configure(state='disabled')

    def download():
        try:
            # Start the 'docker pull' command and stream its output
            docker_process = subprocess.Popen(
                ["docker", "pull", "murphydaviducl/runexomedepth:latest"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1
            )

            # Stream the output
            stream_output(docker_process)

            # Wait for the process to complete and check the return code
            docker_process.wait()
            if docker_process.returncode == 0:
                container_status_label.config(text="Ready to Run", foreground="green")
                status_label.config(text="Container Downloaded")
                global container_downloaded
                container_downloaded = True
            else:
                container_status_label.config(text="Download Failed", foreground="red")
                status_label.config(text="Error")
                console_output.insert(tk.END, "Download failed. See error above.\n")

        except Exception as e:
            container_status_label.config(text="Download Failed", foreground="red")
            status_label.config(text="Error")
            console_output.configure(state='normal')
            console_output.insert(tk.END, f"Error: {str(e)}\n")
            console_output.see(tk.END)
            console_output.configure(state='disabled')

    # Run the download in a separate thread to keep the UI responsive
    threading.Thread(target=download).start()
    
# Initialize the main window
root = tk.Tk()
root.title("ExomeDepth Docker Runner")
root.geometry("700x600")
root.configure(bg='#f0f0f0')  # Light background color for a modern feel

# Style definitions
style = ttk.Style()
style.configure('TButton', font=('Helvetica', 10), padding=5)
style.configure('TLabel', font=('Helvetica', 10), padding=5)
style.configure('TEntry', padding=5)

# Define padding to be used consistently
pad = {'padx': 10, 'pady': 10}

# Create a frame to organize the layout
frame = ttk.Frame(root, padding=20)
frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Test Samples
ttk.Label(frame, text="Test Samples (.bam files):").grid(row=0, column=0, sticky=tk.W, **pad)
test_samples_label = ttk.Label(frame, text="No files selected")
test_samples_label.grid(row=0, column=1, **pad)
ttk.Button(frame, text="Browse", command=select_test_samples).grid(row=0, column=2, **pad)

# Baseline Samples
ttk.Label(frame, text="Baseline Samples (.bam files):").grid(row=1, column=0, sticky=tk.W, **pad)
baseline_samples_label = ttk.Label(frame, text="No files selected")
baseline_samples_label.grid(row=1, column=1, **pad)
ttk.Button(frame, text="Browse", command=select_baseline_samples).grid(row=1, column=2, **pad)

# Target File (optional)
ttk.Label(frame, text="Target File (Optional):").grid(row=2, column=0, sticky=tk.W, **pad)
targets_entry = ttk.Entry(frame, width=50)
targets_entry.grid(row=2, column=1, **pad)
ttk.Button(frame, text="Browse", command=lambda: browse_file(targets_entry)).grid(row=2, column=2, **pad)

# Annotation File (optional)
ttk.Label(frame, text="Annotation File (Optional):").grid(row=3, column=0, sticky=tk.W, **pad)
annotation_entry = ttk.Entry(frame, width=50)
annotation_entry.grid(row=3, column=1, **pad)
ttk.Button(frame, text="Browse", command=lambda: browse_file(annotation_entry)).grid(row=3, column=2, **pad)

# Output Directory
ttk.Label(frame, text="Output Directory:").grid(row=4, column=0, sticky=tk.W, **pad)
output_dir_entry = ttk.Entry(frame, width=50)
output_dir_entry.grid(row=4, column=1, **pad)
ttk.Button(frame, text="Browse", command=lambda: browse_directory(output_dir_entry)).grid(row=4, column=2, **pad)

# Run Button
ttk.Button(frame, text="Run", command=run_docker, style='TButton').grid(row=5, column=1, pady=20)

# Download Container Button
ttk.Button(frame, text="Download Container", command=download_container, style='TButton').grid(row=5, column=2, pady=20)

# Console Output for Logs/Status
console_output_frame = ttk.LabelFrame(root, text="Console Output", padding=(10, 10))
console_output_frame.grid(row=6, column=0, padx=20, pady=10, sticky=tk.W + tk.E + tk.N + tk.S)

console_output = tk.Text(console_output_frame, height=10, width=80, wrap=tk.WORD, bg='#eaeaea', font=('Courier New', 10), state='disabled')
console_output.grid(row=6, column=0, padx=10, pady=10)

# Status Label
status_label = ttk.Label(root, text="Ready")
status_label.grid(row=7, column=0, padx=20, pady=5)

# Docker Status Indicator
docker_status_label = ttk.Label(root, text="Checking Docker...", foreground="orange")
docker_status_label.grid(row=8, column=0, padx=20, pady=5)

# Container Status Indicator
container_status_label = ttk.Label(root, text="Checking Container...", foreground="orange")
container_status_label.grid(row=9, column=0, padx=20, pady=5)

# Handle window close
root.protocol("WM_DELETE_WINDOW", on_closing)

# Check Docker installation and container status
docker_installed = check_docker_installed()
container_downloaded = False
if docker_installed:
    container_downloaded = check_container_downloaded()

root.update_idletasks()
width = root.winfo_reqwidth()
height = root.winfo_reqheight() 
# Set the window size dynamically based on the frame's size
root.geometry(f"{width + 20}x{height + 20}")  # Adding some padding

# Run the Tkinter event loop
root.mainloop()
