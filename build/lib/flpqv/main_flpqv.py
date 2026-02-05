"""
FLPQViewer

Most of parts of ploting figures are managed in plot_figures.py.
This code is written for interactive plot 


"""

import os
import pprint
import copy
import numpy as np

import argparse
import time
import multiprocessing

from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler

import matplotlib.pyplot as plt


from show_figures import show_figures


def read_toml(toml_file):
    
    #  Python バージョンによる `tomllib` or `toml` の切り替え
    try:
        if sys.version_info >= (3, 11):
            import tomllib
            def load_toml(toml_file):
                with open(toml_file, "rb") as f:
                    return tomllib.load(f)
        else:
            import toml
            def load_toml(toml_file):
                with open(toml_file, "r", encoding="utf-8") as f:
                    return toml.load(f)
    except ImportError:
        print("`toml` module cannot be found\n `pip install toml` ")
        sys.exit(1)
    
    #  TOML ファイルを読み込む
    try:
        #print("start read toml_file")
        data_toml = load_toml(toml_file)
        #print("Succeeded reading TOML file")
        #print(data_toml)
        return data_toml
    except FileNotFoundError:
        print(f" cannot find '{toml_file}' ")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    #  TOML ファイルの読み込み終了



def deep_update(original, new):
    for key, value in new.items():
        if isinstance(value, dict) and key in original:
            # originalにそのキーがあって、値が辞書なら再帰的にマージ
            deep_update(original[key], value)
        else:
            # それ以外は上書き
            original[key] = value
    return original

check_update_flag = True
# ファイル変更を検知するクラス
class TomlChangeHandler(FileSystemEventHandler):
    def __init__(self, queue, process2, toml_file, default_config, data_toml):
        super().__init__()
        self.toml_file = toml_file
        self.default_config = default_config
        self.queue = queue
        self.process = process2

    def on_modified(self, event):
        global check_update_flag
        if event.src_path.endswith(self.toml_file):
            if (check_update_flag) :
                user_config = read_toml(self.toml_file)
                data_toml = deep_update(self.default_config, user_config)
                print("toml file is updated")

                self.process.terminate()
                self.process.join()
                self.process = multiprocessing.Process(target=show_figures, args=(data_toml,))
                self.process.start()

                check_update_flag = False
            else:
                check_update_flag = True


# 監視を開始する関数
def start_watching(queue, process2, toml_file, default_config, data_toml):

    event_handler = TomlChangeHandler(queue, process2, toml_file, default_config, data_toml)
    observer = Observer()
    observer.schedule(event_handler, ".", recursive=False)
    observer.start()
    print("Start watching.")
    
    try:
        while True:
            if not process2.is_alive():  # プロセスが終了していたらループを抜ける
                print("Process2 has terminated.")
                observer.stop()
                break

            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()
    process2.terminate()
    process2.join()


#def plot_figures_interactive(queue, process2):
#    while True:
#        try:
#            # キューからデータを取得
#            print("waiting queue_get")
#            data_toml = queue.get()
#            process2.terminate()
#            process2.join()
#            process2 = multiprocessing.Process(target=plot_figures, args=(data_toml,))
#            process2.start()
#            #plot_figures(data_toml)      
#
#            #data_toml = queue.get()
#            #process2.terminate()
#            ##data_toml = queue.get(timeout=1)
#            #process2 = multiprocessing.Process(target=plot_figures_interactive, args=(queue,))
#            #process2.start()
#            ##plt.close()
#            #plot_figures(data_toml)
#            #if (data_toml["show_plot"]):
#            #    plt.show(block=False)
#
#        except queue.Empty:
#            continue  # キューが空ならスキップ

if __name__ == '__main__':

    import sys

    parser = argparse.ArgumentParser(description="FLPQViewer for visualization of FLPQ data.")
    parser.add_argument("filename", nargs="?", help="File to read (optional).")
    parser.add_argument("-cite", action="store_true", help="Print the BibTeX citations.")
    args = parser.parse_args()
    
    if args.cite:
        print(get_bibtex())
        exit()
    
    ##  コマンドライン引数のチェック
    #if len(sys.argv) < 2:
    #    print("Specify TOML file name")
    #    print("How to use: python script.py config.toml")
    #    sys.exit(1)


    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(script_dir + "/default.toml")
    default_config =  read_toml(script_dir + "/default.toml")
    
    #toml_file = sys.argv[1]  # コマンドライン引数からファイル名を取得
    toml_file = args.filename  # コマンドライン引数からファイル名を取得


    # 初回表示
    user_config = read_toml(toml_file)
    #print(data_toml)

    data_toml = deep_update(default_config, user_config)


    if(data_toml["save"]["interactive"]):
        

        queue = multiprocessing.Queue()

        multiprocessing.set_start_method("spawn", force=True)

        process2 = multiprocessing.Process(target=show_figures, args=(data_toml,))
        process2.start()

        start_watching(queue, process2, toml_file, default_config, data_toml)


    try:
        if( data_toml["save"]["pdf"] or
            data_toml["save"]["png"] or
            data_toml["save"]["Plotly"] or
            data_toml["save"]["obj"]):

            data_toml["save"]["interactive"] = False
            show_figures(data_toml)
    except KeyError:
        pass

    sys.exit(0)

    #plot_isosurface("LPQ.kinetic_energy_density.cube", "LPQ.chemical_potential.cube", max_value=-0.20, min_value=-0.38, isovalue=0.00001, contour_levels=5)
    #plot_isosurface("LPQ.kinetic_energy_density.cube", "LPQ.chemical_potential.cube", max_value=-0.175, min_value=-0.2, isovalue=0.00001, contour_levels=5)
    #plot_isosurface("RESULT.tden.cube", "RESULT.vhart.cube", max_value=0.02, min_value=-0.02, isovalue=0.0005)
