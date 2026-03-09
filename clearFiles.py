import os

folder_paths = ['/home/peptool/mysite/static/serve/images', '/home/peptool/mysite/static/serve/files']

for folder_path in folder_paths:
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

print("Files cleared successfully")
