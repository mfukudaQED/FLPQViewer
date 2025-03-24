import re

Mypackage = "flpqv"

def get_current_version(file_path):
    """指定したファイルから現在のバージョンを取得する関数"""
    with open(file_path, 'r') as file:
        content = file.read()

    # 正規表現でバージョンを抽出（例：__version__ = "1.0.0"）
    match = re.search(r'__version__ = "([^"]+)"', content)
    if match:
        return match.group(1)
    else:
        raise ValueError("Version not found in the file!")

def update_version(file_path, new_version):
    """指定したファイルでバージョンを更新する関数"""
    with open(file_path, 'r') as file:
        content = file.read()

    # バージョンを正規表現で探して置き換え
    content = re.sub(r'__version__ = "[^"]+"', f'__version__ = "{new_version}"', content)

    with open(file_path, 'w') as file:
        file.write(content)

def update_pyproject_version(file_path, new_version):
    """pyproject.toml内のバージョンを更新する関数"""
    with open(file_path, 'r') as file:
        content = file.read()

    # pyproject.toml内のバージョンを置き換え
    content = re.sub(r'version = "[^"]+"', f'version = "{new_version}"', content)

    with open(file_path, 'w') as file:
        file.write(content)

if __name__ == "__main__":
    # 現在のバージョンを取得
    current_version = get_current_version(f"{Mypackage}/__init__.py")
    print(f"Current version: {current_version}")

    # ユーザーに新しいバージョンを入力させる
    new_version = input("Enter the new version: ")

    # 新しいバージョンを更新
    update_version(f"{Mypackage}/__init__.py", new_version)
    update_pyproject_version("pyproject.toml", new_version)

    print(f"Version updated from {current_version} to {new_version}")
