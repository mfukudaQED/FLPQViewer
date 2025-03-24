from mayavi import mlab

mlab.test_plot3d()  # 何かプロットする
mlab.show()  # ← これがないとウィンドウがすぐ閉じる！


#import numpy as np
#from mayavi import mlab
#
## グリッドの作成
#x = np.linspace(-1, 1, 20)
#y = np.linspace(-1, 1, 20)
#X, Y = np.meshgrid(x, y)
#
## 1つ目の平面: z = x + y
#Z1 = X + Y
#
## 2つ目の平面: z = 1 - x
#Z2 = 1 - X
#
## Mayaviで可視化
#mlab.figure(bgcolor=(1, 1, 1))  # 背景色を白に
#
## 1つ目の平面
#mlab.mesh(X, Y, Z1, color=(1, 0, 0), opacity=0.6)  # 赤色
#
## 2つ目の平面
#mlab.mesh(X, Y, Z2, color=(0, 0, 1), opacity=0.6)  # 青色
#
## 軸の表示
#mlab.axes(xlabel='X', ylabel='Y', zlabel='Z')
#
## 表示
#mlab.show()
