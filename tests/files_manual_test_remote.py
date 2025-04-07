#!/usr/bin/env python3
"""
File Service API Integration Test - run manually with python

This script tests the DeepOrigin FilesClient API by performing a series of
operations against a remote file service:

1. Creates random test files locally in temp_test_data/
2. Tests upload_file functionality for individual files
3. Tests get_metadata to verify file properties
4. Tests download_file to retrieve files
5. Tests delete_file functionality
6. Tests sync_dir for bidirectional synchronization
7. Cleans up by removing remote test files

If the script is run and temp_test_data already exists, it will be deleted
and recreated to ensure a clean test environment.
"""

import os
import sys
import shutil
import random
import string
import time
import hashlib
from pathlib import Path
import tempfile
from typing import List, Dict, Any

from deeporigin.files import FilesClient

# ==== Configuration ====
# Define the Base URL or AUTH_TOKEN to override default URL and token the file service API
# If not defined, these come from 'deeporigin config' settings
#API_BASE_URL = "https://os.edge.deeporigin.io"
#AUTH_TOKEN = "YOUR_AUTH_TOKEN"

# Test directories
TEST_ROOT_DIR = "temp_test_data"
SYNC_SOURCE_DIR = os.path.join(TEST_ROOT_DIR, "sync_data")
SYNC_SUBFOLDER = os.path.join(SYNC_SOURCE_DIR, "subfolder")
SYNC_DOWNLOAD_DIR = os.path.join(TEST_ROOT_DIR, "sync_download")

# Remote path prefix (where files will be stored in remote service)
REMOTE_PREFIX = "files:///test_files"
FILE_PREFIX = "td_"  # Prefix for test data files

# Number of files to create
NUM_SYNC_FILES = 10  # Files in sync_data dir
NUM_SUBFOLDER_FILES = 2  # Files in subfolder
NUM_UPLOAD_TEST_FILES = 2  # Files to explicitly test upload/download

# File size range in bytes
MIN_FILE_SIZE = 2 * 1024  # 2KB
MAX_FILE_SIZE = 30 * 1024  # 30KB

# ==== Helper functions ====
def generate_random_text(size: int) -> str:
    """
    Generate random text content with specified size, wrapped at 80 chars.    
    Args:
        size: Desired size in bytes        
    Returns:
        Random string of specified size
    """
    chars = string.ascii_letters + string.digits + string.punctuation + ' ' * 10
    lines = []
    current_line = []
    total_size = 0
    
    while total_size < size:
        if len(''.join(current_line)) >= 80:
            lines.append(''.join(current_line))
            current_line = []
        
        char = random.choice(chars)
        current_line.append(char)
        total_size += len(char.encode('utf-8'))
    
    if current_line:
        lines.append(''.join(current_line))
    
    return '\n'.join(lines)

def create_test_file(path: str, size: int = None) -> int:
    """
    Create a test file with random content.    
    Args:
        path: Path where file should be created
        size: Size of file to create (random if None)        
    Returns:
        Size of created file
    """
    if size is None:
        size = random.randint(MIN_FILE_SIZE, MAX_FILE_SIZE)
        
    content = generate_random_text(size)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    # Return actual size
    return os.path.getsize(path)

def get_file_hash(file_path: str) -> str:
    """
    Calculate SHA-256 hash of a file.    
    Args:
        file_path: Path to file        
    Returns:
        Hex digest of file hash
    """
    hash_obj = hashlib.sha256()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_obj.update(chunk)
    return hash_obj.hexdigest()

def print_header(text: str):
    """Print a formatted header for test sections."""
    print("\n" + "=" * 80)
    print(f"  {text}")
    print("=" * 80)

def print_step(text: str):
    """Print a formatted step message."""
    print(f"\n--> {text}")

def print_remote_files(client, remote_path):
    """
    List and print files at the given remote path.    
    Args:
        client: FilesClient instance
        remote_path: Remote path to list files from
    """
    print_step(f"Listing files at {remote_path}")
    
    try:
        files = client.list_dir(remote_path)
        if files and len(files) > 0:
            print(f"  Found {len(files)} files/directories:")
            for file_info in files:
                # Extract relevant information from file_info
                # Adjust these fields based on the actual structure returned by list_dir
                name = file_info.get("Key", "Unknown")
                size = file_info.get("Size", "Unknown")
                #last_modified = file_info.get("LastModified", "Unknown")
                file_type = "Directory" if file_info.get("is_dir", False) else "File"
                
                print(f"  - {name} ({file_type}, {size} bytes)")
        else:
            print("  No files found or empty directory")
    except Exception as e:
        print(f"  Failed to list files: {e}")


# ==== Main test function ====
def run_tests():
    """Run the file service integration tests."""
    print_header("File Service API Integration Test")

    # Initialize client    
    try:
        client = FilesClient(
            base_url=API_BASE_URL if 'API_BASE_URL' in globals() else None,
            token=AUTH_TOKEN if 'AUTH_TOKEN' in globals() and AUTH_TOKEN else None
        )
        
        if not client.check_health():
            print(f"❌ File service health check failed")
            print(f"Please check if the service at {API_BASE_URL} is running")
            sys.exit(1)        
        if not client.check_auth_ok():
            print(f"❌ Authentication failed")
            print(f"Please check your AUTH_TOKEN")
            sys.exit(1)
            
        print("✅ Authentication successful")
        
    except Exception as e:
        print(f"❌ Failed to establish connection: {str(e)}")
        print(f"Please check your API_BASE_URL ({API_BASE_URL}) and AUTH_TOKEN")
        sys.exit(1)  # Exit with error code

    # Clean up any existing test data
    if os.path.exists(TEST_ROOT_DIR):
        print_step(f"Removing existing {TEST_ROOT_DIR} directory")
        shutil.rmtree(TEST_ROOT_DIR)
    
    # Create test directories
    print_step("Creating test directories")
    os.makedirs(SYNC_SUBFOLDER, exist_ok=True)
    os.makedirs(SYNC_DOWNLOAD_DIR, exist_ok=True)    
    
    # Create test files for explicit upload/download testing
    print_step(f"Creating {NUM_UPLOAD_TEST_FILES} test files for upload/download tests")
    upload_test_files = []
    for i in range(NUM_UPLOAD_TEST_FILES):
        file_path = os.path.join(TEST_ROOT_DIR, f"{FILE_PREFIX}upload_test_{i}.txt")
        file_size = create_test_file(file_path)
        upload_test_files.append({"path": file_path, "size": file_size})
        print(f"  Created {file_path} ({file_size} bytes)")
    
    # Create test files for sync testing
    print_step(f"Creating {NUM_SYNC_FILES} test files in sync directory")
    sync_files = []
    for i in range(NUM_SYNC_FILES):
        file_path = os.path.join(SYNC_SOURCE_DIR, f"{FILE_PREFIX}sync_{i}.txt")
        file_size = create_test_file(file_path)
        sync_files.append({"path": file_path, "size": file_size})
        print(f"  Created {file_path} ({file_size} bytes)")
    
    print_step(f"Creating {NUM_SUBFOLDER_FILES} test files in subfolder")
    subfolder_files = []
    for i in range(NUM_SUBFOLDER_FILES):
        file_path = os.path.join(SYNC_SUBFOLDER, f"{FILE_PREFIX}subfolder_{i}.txt")
        file_size = create_test_file(file_path)
        subfolder_files.append({"path": file_path, "size": file_size})
        print(f"  Created {file_path} ({file_size} bytes)")
    
    # Define remotes paths for test files
    all_files = upload_test_files + sync_files + subfolder_files
    for file_info in all_files:
        local_path = file_info["path"]
        relative_path = os.path.relpath(local_path, TEST_ROOT_DIR)
        file_info["remote_path"] = f"{REMOTE_PREFIX}/{relative_path}"
    

    # === Test 1: Upload files ===
    print_header("Test 1: Upload files")
    upload_success = True
    
    for file_info in upload_test_files:
        local_path = file_info["path"]
        remote_path = file_info["remote_path"]
        print_step(f"Uploading {local_path} to {remote_path}")

        try:
            success = client.upload_file(src=local_path, dest=remote_path)
            if success:
                print(f"  Upload successful")
                file_info["uploaded"] = True
            else:
                print(f"  Upload failed")
                upload_success = False
                file_info["uploaded"] = False
        except Exception as e:
            print(f"  Upload failed with error: {e}")
            upload_success = False
            file_info["uploaded"] = False


    # === Test 2: Get metadata ===
    print_header("Test 2: Get metadata")
    
    if not upload_success:
        print("Skipping metadata test as upload failed")
    else:
        for file_info in upload_test_files:
            if not file_info.get("uploaded", False):
                continue
                
            remote_path = file_info["remote_path"]
            expected_size = file_info["size"]
            print_step(f"Getting metadata for {remote_path}")
            
            try:
                metadata = client.get_metadata(remote_path)
                if metadata:
                    print(f"  Metadata retrieved: {metadata}")
                    
                    # Try to extract file size from metadata
                    # This depends on the actual API response format
                    file_size = None
                    if "content-length" in metadata:
                        file_size = int(metadata["content-length"])
                    elif "Content-Length" in metadata:
                        file_size = int(metadata["Content-Length"])
                        
                    if file_size is not None:
                        if file_size == expected_size:
                            print(f"  Size matches: {file_size} bytes")
                            file_info["size_verified"] = True
                        else:
                            print(f"  Size mismatch: expected {expected_size}, got {file_size}")
                            file_info["size_verified"] = False
                    else:
                        print("  Could not find file size in metadata")
                        file_info["size_verified"] = False
                else:
                    print("  Failed to retrieve metadata")
                    file_info["size_verified"] = False
            except Exception as e:
                print(f"  Metadata retrieval failed with error: {e}")
                file_info["size_verified"] = False
    
    
    # === Test 3: Download files ===
    print_header("Test 3: Download files")    

    if not upload_success:
        print("Skipping download test as upload failed")
    else:
        for file_info in upload_test_files:
            if not file_info.get("uploaded", False):
                continue
                
            remote_path = file_info["remote_path"]
            original_path = file_info["path"]
            
            # Create a new filename for the downloaded version
            base_name = os.path.basename(original_path)
            name_parts = os.path.splitext(base_name)
            download_filename = f"{name_parts[0]}_download{name_parts[1]}"
            download_path = os.path.join(TEST_ROOT_DIR, download_filename)
            
            print_step(f"Downloading {remote_path} to {download_path}")
            
            try:
                success = client.download_file(src=remote_path, dest=download_path)
                if success:
                    print(f"  Download successful")
                    
                    # Verify content integrity
                    original_hash = get_file_hash(original_path)
                    download_hash = get_file_hash(download_path)
                    
                    if original_hash == download_hash:
                        print("  Content verification passed: files are identical")
                        file_info["download_verified"] = True
                    else:
                        print("  Content verification failed: files are different")
                        file_info["download_verified"] = False
                else:
                    print(f"  Download failed")
                    file_info["download_verified"] = False
            except Exception as e:
                print(f"  Download failed with error: {e}")
                file_info["download_verified"] = False
    

    # === Test 4: Delete files ===
    print_header("Test 4: Delete files")
    
    if not upload_success or len(upload_test_files) < 1:
        print("Skipping delete test as upload failed or no files available")
    else:
        # Delete one of the test files
        file_to_delete = upload_test_files[0]
        remote_path = file_to_delete["remote_path"]
        
        print_step(f"Deleting {remote_path}")
        
        try:
            success = client.delete_file(remote_path)
            if success:
                print(f"  Delete successful")
                file_to_delete["deleted"] = True
                
                # Verify deletion by trying to get metadata
                print_step(f"Verifying deletion by checking metadata")
                try:
                    metadata = client.get_metadata(remote_path)
                    if metadata:
                        print("  File still exists - deletion verification failed")
                        file_to_delete["deletion_verified"] = False
                    else:
                        print("  File not found - deletion verification successful")
                        file_to_delete["deletion_verified"] = True
                except Exception as e:
                    # If we get an error (like 404), the file was likely deleted
                    print(f"  Got error when checking metadata (expected): {e}")
                    file_to_delete["deletion_verified"] = True
            else:
                print(f"  Delete failed")
                file_to_delete["deleted"] = False
        except Exception as e:
            print(f"  Delete failed with error: {e}")
            file_to_delete["deleted"] = False
    

    # === Test 5: Test sync_dir (upload) ===
    print_header("Test 5: Test sync_dir (upload)")
    
    sync_remote_path = f"{REMOTE_PREFIX}/sync_data"
    print_step(f"Syncing {SYNC_SOURCE_DIR} to {sync_remote_path}")
    
    try:
        success = client.sync_dir(src=SYNC_SOURCE_DIR, dst=sync_remote_path)
        if success:
            print(f"  Sync upload successful")
            sync_upload_success = True
            
            # List remote files after successful sync
            print_remote_files(client, sync_remote_path)
            
            # Also list the subfolder
            subfolder_remote_path = f"{sync_remote_path}/subfolder"
            print_remote_files(client, subfolder_remote_path)
        else:
            print(f"  Sync upload failed")
            sync_upload_success = False
    except Exception as e:
        print(f"  Sync upload failed with error: {e}")
        sync_upload_success = False
    

    # === Test 6: Test sync_dir (download) ===
    print_header("Test 6: Test sync_dir (download)")

    # Initialize verification_success at the beginning of the test block
    verification_success = False

    if not sync_upload_success:
        print("Skipping sync download test as sync upload failed")
    else:
        print_step(f"Syncing {sync_remote_path} to {SYNC_DOWNLOAD_DIR}")
                
        try:
            success = client.sync_dir(src=sync_remote_path, dst=SYNC_DOWNLOAD_DIR)
            if success:
                print(f"  Sync download successful")
                
                # Verify content integrity for a sample of files
                print_step("Verifying content integrity of synced files")

                verification_success = True
                
                for file_info in sync_files:  # Check files
                    local_path = file_info["path"]
                    relative_path = os.path.relpath(local_path, SYNC_SOURCE_DIR)
                    download_path = os.path.join(SYNC_DOWNLOAD_DIR, relative_path)
                    
                    if os.path.exists(download_path):
                        original_hash = get_file_hash(local_path)
                        download_hash = get_file_hash(download_path)
                        
                        if original_hash == download_hash:
                            print(f"  Content verification passed for {relative_path}")
                        else:
                            print(f"  Content verification failed for {relative_path}")
                            verification_success = False
                    else:
                        print(f"  File not found: {download_path}")
                        verification_success = False
                
                if verification_success:
                    print("  All checked files verified successfully")
                else:
                    print("  Some files failed verification")
            else:
                print(f"  Sync download failed")
        except Exception as e:
            print(f"  Sync download failed with error: {e}")
    

    # === Test 7: Cleanup - Delete all remaining remote files ===
    print_header("Test 7: Cleanup")
    
    # Only attempt to delete files that were successfully uploaded
    files_to_cleanup = [file_info for file_info in all_files 
                        if file_info.get("uploaded", False) and not file_info.get("deleted", False)]
    
    if not files_to_cleanup and not sync_upload_success:
        print("Skipping cleanup as no files were successfully uploaded")
    else:
        for file_info in files_to_cleanup:
            remote_path = file_info["remote_path"]
            print_step(f"Deleting {remote_path}")
            
            try:
                success = client.delete_file(remote_path)
                if success:
                    print(f"  Delete successful")
                else:
                    print(f"  Delete failed")
            except Exception as e:
                print(f"  Delete failed with error: {e}")
        
        # Only attempt to delete sync directory if sync was successful
        if sync_upload_success:
            sync_dir_path = f"{REMOTE_PREFIX}/sync_data"
            print_step(f"Attempting to delete directory: {sync_dir_path}")
            try:
                success = client.delete_file(sync_dir_path)
                if success:
                    print(f"  Directory delete successful")
                else:
                    print(f"  Directory delete failed")
            except Exception as e:
                print(f"  Directory delete failed with error: {e}")
        
    
    # Initialize all test results as failed by default
    upload_test_result = "FAILED"
    metadata_test_result = "FAILED"
    download_test_result = "FAILED"
    delete_test_result = "FAILED"
    sync_upload_test_result = "FAILED"
    sync_download_test_result = "FAILED"
    
    # Only mark as PASSED if we have evidence of success
    if upload_success:
        upload_test_result = "PASSED"
        
        # Only check these if upload succeeded
        if all(file_info.get("size_verified", False) for file_info in upload_test_files if file_info.get("uploaded", False)):
            metadata_test_result = "PASSED"
        
        if all(file_info.get("download_verified", False) for file_info in upload_test_files if file_info.get("uploaded", False)):
            download_test_result = "PASSED"
        
        if any(file_info.get("deletion_verified", False) for file_info in upload_test_files):
            delete_test_result = "PASSED"
    
    if sync_upload_success:
        sync_upload_test_result = "PASSED"
        
        # Only check this if sync upload succeeded
        if verification_success:
            sync_download_test_result = "PASSED"


    # Clean up any existing test data
    if (upload_test_result == "PASSED" and metadata_test_result == "PASSED" and download_test_result == "PASSED" \
    and delete_test_result == "PASSED" and sync_upload_test_result == "PASSED" and sync_download_test_result == "PASSED"):
        if os.path.exists(TEST_ROOT_DIR):
            print_step(f"Removing existing {TEST_ROOT_DIR} directory")
            shutil.rmtree(TEST_ROOT_DIR)
    else:
        print("\nTest data is available in:", TEST_ROOT_DIR)    

    print_header("Test Summary")    
    
    print("1. Upload test:", upload_test_result)
    print("2. Metadata test:", metadata_test_result)
    print("3. Download test:", download_test_result)
    print("4. Delete test:", delete_test_result)
    print("5. Sync upload test:", sync_upload_test_result)
    print("6. Sync download test:", sync_download_test_result)
    
if __name__ == "__main__":
    run_tests() 