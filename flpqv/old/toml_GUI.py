import argparse
import tkinter as tk
from tkinter import messagebox
import toml

# TOMLファイルの読み込み
def load_toml(file_path):
    try:
        with open(file_path, 'r') as file:
            return toml.load(file)
    except Exception as e:
        messagebox.showerror("Error", f"Error reading the TOML file: {e}")
        return {}

# TOMLファイルの保存
def save_toml(file_path, data):
    try:
        with open(file_path, 'w') as file:
            toml.dump(data, file)
        messagebox.showinfo("Success", "File saved successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"Error saving the TOML file: {e}")

# GUIの設定
def create_gui(toml_data, file_path):
    def update_toml():
        # 入力された内容をTOMLデータに反映
        for section, values in toml_data.items():
            if isinstance(values, dict):  # valuesが辞書型である場合に処理
                for key, entry in entries[section].items():
                    if isinstance(entry, tk.StringVar):
                        toml_data[section][key] = entry.get()  # 通常のエントリ
                    elif isinstance(entry, tk.BooleanVar):
                        toml_data[section][key] = entry.get()  # チェックボックスの値
        save_toml(file_path, toml_data)

    window = tk.Tk()
    window.title("TOML Editor")
    window.configure(bg="#f5f5f5")

    # スクロール用のフレーム
    canvas = tk.Canvas(window)
    canvas.pack(side="left", fill="both", expand=True)

    scrollbar = tk.Scrollbar(window, orient="vertical", command=canvas.yview)
    scrollbar.pack(side="right", fill="y")

    canvas.configure(yscrollcommand=scrollbar.set)
    scrollable_frame = tk.Frame(canvas, bg="#f5f5f5")

    # スクロール可能なフレームを作成
    scrollable_frame.pack(fill="both", expand=True)

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")


    entries = {}
    row = 0

    # セクションごとにウィジェットを作成
    for section, section_data in toml_data.items():
        if isinstance(section_data, dict):  # セクションデータが辞書型の場合
            # セクション名を表示
            label = tk.Label(scrollable_frame, text=f"[{section}]", bg="#f5f5f5", anchor="w", padx=10, font=("Arial", 12, "bold"))
            label.pack(anchor="w", padx=10, pady=5, fill='x')
            row += 1

            section_entries = {}
            for key, value in section_data.items():
                label = tk.Label(scrollable_frame, text=key, bg="#f5f5f5", anchor="w", padx=10, font=("Arial", 12))
                label.pack(anchor="w", padx=10, pady=5, fill='x')

                if isinstance(value, bool):  # True/False の場合はチェックボックス
                    var = tk.BooleanVar(value=value)
                    entry = tk.Checkbutton(scrollable_frame, variable=var, bg="#f5f5f5", activebackground="#f5f5f5")
                    entry.pack(anchor="w", padx=10, pady=5, fill='x')
                    section_entries[key] = var
                else:  # それ以外の値はエントリ
                    entry = tk.Entry(scrollable_frame, relief="flat", bd=0, highlightthickness=1, highlightcolor="#3498db", highlightbackground="#ccc", width=40)
                    entry.pack(anchor="w", padx=10, pady=5, fill='x')
                    entry.insert(0, str(value))
                    section_entries[key] = entry

                row += 1

            entries[section] = section_entries

    # 保存ボタン
    save_button = tk.Button(scrollable_frame, text="Save", command=update_toml, bg="#3498db", fg="white", relief="flat", padx=15, pady=8)
    save_button.pack(pady=10)

    window.mainloop()



if __name__ == '__main__':

    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", nargs="?", help="File to read (optional).")
    args = parser.parse_args()


    # TOMLファイルを読み込んでGUIを表示
    file_path = args.filename
    toml_data = load_toml(file_path)
    create_gui(toml_data, file_path)
