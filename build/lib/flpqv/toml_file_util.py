"""
TOML設定ファイルユーティリティ
"""
from logging import exception
from typing import Any
import tkinter as tk
import sys
import toml
import tkinter.font as tkfont

fontname = "gothic"
fontsize = 11

class TomlFileUtil:
    """
    tomlファイル操作、設定画面サポートユーティリティ
    """
    def __init__(self) -> None:
        self.path = ""      # tomlファイルのパス
        self.toml_doc = {}  # tomlファイルを読んで作成した辞書
        self._var_dict = {} # Tkinter variable変数(ウィジェット変数)を値にした辞書

    def dump_toml(self, toml_dict:dict, path:str):
        """
        tomlで読み込んだ辞書をtomlファイルに出力する
		Args:
            dict:   tomlで読み込んだ辞書
            str:    保存先ファイル名
        """
        try:
            with open(path, 'w', encoding="utf-8") as configfile:
                toml.dump(toml_dict, configfile)
        except exception as e:
            print(e)

    def save_toml(self, save_path:str="", section:str="", event=None):
        """
        ウィジェット変数の値をsectionで指定したセクションに移してsave_pathで指定したファイルに保存
		Args:
            str:    保存先ファイル名
            str:    保存先セクション
        """
        if not save_path:
            return
        self.set_var_dict2toml_dict(section)
        self.dump_toml(self.toml_doc, save_path)

    def load_default(self, event=None):
        """
        ウィジェット変数の辞書にTOMLのDEFAULTセクションの値を設定
        """
        self.set_toml2var_dict("DEFAULT")

    def read_toml(self, path:str) -> dict:
        """
        tomlファイルを読み込み辞書を返す
		Args:
            str:    読み込むファイル名
		Returns:
			dict:   tomlファイルの辞書。エラーの時はNone
        """
        self.path = path
        try:
            with open(path, "r", encoding="utf-8") as f:
                toml_dict = toml.load(f)
        except Exception as e:
            print(e)
            return None
        self.toml_doc = toml_dict
        return toml_dict

    def entry_validate(self, action:str, modify_str:str) -> bool:
        """
        エントリーの入力検証
        Args:
            str:    アクション(削除：0、挿入：1、その他：-1)
            str:    挿入、削除されるテキスト
        """
        if action != "1": return True   # 挿入の時だけ検証
        return modify_str.isdigit()     # 数字かどうか

    def create_frame_from_toml_dict(self, parent:Any
        , has_save_button:bool=False, has_default_button:bool=False, is_grid:bool=False) -> dict:
        """
        tomlデータからDEFINITIONセクション(キー:[説明, 型]の辞書)を抜き出し、フレームを作成
        tomlデータは事前に読んでおく
        保存ボタンを表示する場合、ボタンに割り付けるメソッドを後から指定する
            例 obj.btn_save.config(command=lambda path=my_path: toml.save_toml(path, "USER"))
		Args:
            Any:    主にFrame 親コンテナ
            bool:   保存用ボタンの有無
            bool:   labelとentryを横並びにするかどうか
		Returns:
			dict:   ウィジェット変数の辞書
        """
        toml_dict = self.toml_doc["DEFINITION"]
        self._var_dict = {}
        # 説明の最大文字数を求める
        max_length = max([len(v2[0]) for v1 in toml_dict.values() for v2 in v1.values()])
        for group_key, group_item in toml_dict.items():
            label_frame = tk.LabelFrame(parent, text=group_key, padx=5, pady=5, font=(fontname, fontsize))
            label_frame.pack(fill="x", padx=5,pady=5)
            for i, (k, item) in enumerate(group_item.items()):
                if item[1] == "str":
                    self._var_dict[k] = tk.StringVar()
                    entry_ = tk.Entry(label_frame, textvariable=self._var_dict[k], font=(fontname, fontsize))
                if item[1] == "int":
                    self._var_dict[k] = tk.IntVar()
                    entry_ = tk.Entry(label_frame, textvariable=self._var_dict[k]
                            , validate="key", vcmd=(parent.register(self.entry_validate), "%d", "%S"), font=(fontname, fontsize))
                if item[1] == "bool":
                    self._var_dict[k] = tk.BooleanVar()
                    checkbutton_ = tk.Checkbutton(label_frame, variable=self._var_dict[k]
                        , text=item[0], anchor="w", width=max_length*2, font=(fontname, fontsize))
                    if is_grid:
                        checkbutton_.grid(row=i, column=0, columnspan=2 , sticky="w")
                    else:
                        checkbutton_.pack()
                    continue
                label_ = tk.Label(label_frame, text=item[0], font=(fontname,fontsize))
                # pack,grid
                if is_grid:
                    label_.grid(row=i, column=0, sticky="e")
                    entry_.grid(row=i, column=1)
                else:
                    label_.pack(anchor="w")
                    entry_.pack()
            if is_grid:
                label_frame.columnconfigure(0, weight=1)    # grid column 0 の余白を引き延ばす
        # 保存ボタンの作成
        if has_save_button:
            self.btn_save = tk.Button(parent, text="設定保存", font=(fontname,fontsize))
            self.btn_save.pack(side="bottom", fill="x")

        # デフォルトロードボタンの作成
        if has_default_button:
            self.btn_default = tk.Button(parent, text="デフォルトに戻す", font=(fontname,fontsize)
                , command=self.load_default)
            self.btn_default.pack(side="bottom", fill="x")

        return self._var_dict

    def set_toml2var_dict(self, section:str):
        """
        tomlの辞書からsectionで指定したセクションをウィジェット変数の辞書に設定する
		Args:
            str:    セクション
        """
        for group_item in self.toml_doc[section].values():
            for k, item in group_item.items():
                self._var_dict[k].set(item)

    def set_var_dict2toml_dict(self, section:str):
        """
        ウィジェット変数の辞書の設定値をtomlの辞書sectionで指定したセクションに設定する
		Args:
            str:    セクション
        """
        for group_item in self.toml_doc[section].values():
            for k in group_item:
                group_item[k] = self._var_dict.get(k).get()

if __name__ == '__main__':

#    root = tk.Tk()
#    
#    # システムにインストールされているフォントファミリーを取得
#    available_fonts = tkfont.families()
#    
#    # 使用したいフォント名
#    font_to_check = "Arial"
#    
#    if font_to_check in available_fonts:
#        print(f"'{font_to_check}' フォントはシステムにインストールされています！")
#    else:
#        print(f"'{font_to_check}' フォントはシステムにインストールされていません。")
#
#    print(tkfont.families())
#    
#    root.mainloop()


    data_toml = TomlFileUtil()

    # ファイルから編集画面を作成し起動
    if len(sys.argv) < 2:  
        print("引数にファイルを指定してください。")
        sys.exit()
    my_path = sys.argv[1]
    toml_doc = data_toml.read_toml(my_path)
    if not toml_doc:
        sys.exit()
    root = tk.Tk()
    var_dict = data_toml.create_frame_from_toml_dict(root, True, True, True)
    data_toml.set_toml2var_dict("USER")
    data_toml.btn_save.config(command=lambda path=my_path: toml.save_toml(path, "USER"))
    root.mainloop()
