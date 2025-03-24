import os
import pprint
import copy
import argparse
from collections import ChainMap
import time
import threading
import queue

import numpy as np
import trimesh

from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.distance import cdist

from skimage.measure import marching_cubes

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler

#
#from parameters import BohrR
BohrR = 0.529177249 # Angstrom

class Plot_style_1d:
    def __init__(self) -> None:
        self.plt_style()
    
    def plt_style(self):
        plt.rcParams['font.family'] ='Times New Roman'
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.0
        plt.rcParams['errorbar.capsize'] = 6
        plt.rcParams['lines.markersize'] = 7
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['mathtext.fontset'] = 'cm'
        self.line_styles = ['-', '--', '-.', ':']
        self.markers = ['o', ',', '.', 'v', '^', '<', '>', '1', '2', '3', '.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3']



class Data_cube:
    """ class of data from md file """
    def __init__(self,
                 atomnum = None,
                 numbers = None,   # int
                 valence = None,   # int
                 positions_bohr = None,  # bohr, coordinate of atoms
                 positions_ang = None,  # Ang, coordinate of atoms
                 origin = None, # Ang
                 num_mesh = None,
                 mesh_vec = None, # Ang primitive vector for mesh
                 mesh_grid = None, # Ang
                 val = None,
                 small_mesh_grid = None, # Ang
                 small_val = None):

        pass

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

    def read(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            #self.origin = [float(coord) for coord in line[1:4]]
            self.origin = [float(coord)*BohrR for coord in line[1:4]]


            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]


            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            for _ in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])

                          ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

            self.val = np.empty((self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2])).tolist()

            # Read the file using the next function until the end of the file
            idx = 0
            while True:
                try:
                    line = next(file).strip().split()
                    for value in line:
                        self.val[idx] = float(value)
                        #print(self.val[idx])
                        idx += 1

                except StopIteration:
                    # End of file reached, break out of the loop
                    break

            ### Data structure ###
            #idx = 0
            #for i in range(self.num_mesh[0]):
            #    x = self.origin[0] + i*self.mesh_vec[0]
            #    for j in range(self.num_mesh[1]):
            #        y = self.origin[1] + j*self.mesh_vec[1]
            #        for k in range(num_mesh[2]):
            #            z = self.origin[2] + k*self.mesh_vec[2]
            #            print(x,y,z,self.val[idx])
            #            idx += 1
        return

    def read_header(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            self.origin = [float(coord)*BohrR for coord in line[1:4]]

            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]

            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            for _ in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])

            ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

        return

class Data_atoms():
    """ class of color list of atoms """
    def __init__(self,
                 name = None,
                 color = None,
                 ):

        self.name =\
        ["E", "H",  "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F",\
        "Ne", "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar",  "K",\
        "Ca", "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",\
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y",\
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",\
        "Sn", "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",\
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",\
        "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au",\
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",\
        "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",\
        "Fm", "Md", "No", "Lr"]

        # number = List_element.index("H")                                                                                                                                                                                                                                                                                         
        self.color =\
        ["#FFFFFF", "#FAF0E6", "#F0E68C", "#EE82EE", "#DB7093", "#56256e", "#778899", "#0033ff", "#ee0011", "#7fff7f",\
         "#FF00FF", "#48D1CC", "#8B008B", "#D8BFD8", "#00CED1", "#BC043C", "#FFFF55", "#66CDAA", "#FF69B4", "#ADD8E6",\
         "#5c9291", "#483D8B", "#8B708B", "#008080", "#B22222", "#a033a0", "#BC8F8F", "#50C878", "#66CDAA", "#D2691E",\
         "#BC8F8F", "#F08080", "#E6E6FA", "#DAA520", "#e0ebaf", "#b44c97", "#00a497", "#b77b57", "#006e54", "#d69090",\
         "#634950", "#44617b", "#E9967A", "#FF7F50", "#FFA07A", "#A0522D", "#FFE4C4", "#C0C0C0", "#BDB76B", "#6B8E23",\
         "#556B2F", "#ADFF2F", "#008000", "#7a4171", "#7FFFD4", "#008080", "#00aaaa", "#999900", "#4682B4", "#4169E1",\
         "#7B68EE", "#8A2BE2", "#BA55D3", "#800080", "#FF1493", "#c89932", "#5c9291", "#9d5b8b", "#eeeaec", "#e29676",\
         "#5f6527", "#c70067", "#e9dacb", "#28ff93", "#0582ff", "#baff75", "#43676b", "#47585c", "#fde8d0", "#FFD700",\
         "#dcd6d9", "#bf794e", "#f5b1aa", "#cd5e3c", "#95859c", "#71686c", "#203744", "#ec6d71", "#b55233", "#a19361",\
         "#cc3399", "#3399cc", "#339966", "#ffff00", "#ccff33", "#cc3300", "#cc6600", "#ff0033", "#660066", "#006666",\
         "#ffcc00", "#33cccc", "#99ffff", "#ff6666"]

        self.size =\
        [0.20, 0.41, 0.41, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
         0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03,
         1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32,
         1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90,
         1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.42, 1.42,
         1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03,
         2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90,
         1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
         1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15,
         2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60]

    def hex_to_rgb(hex_code):
        hex_code = hex_code.lstrip("#")  # "#" を削除
        return tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))
        # print(hex_to_rgb("#FF5733"))  # (255, 87, 51)

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

def plot_structure(ax, data_toml, data_cube):
    """
    原子構造をプロットするコード
    """
    # 原子座標を取得
    positions = data_cube.positions_ang
    numbers = np.int64(data_cube.numbers)

    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
        # **重心を求める**
        centroid = np.mean(data_cube.positions_ang, axis=0)

    positions = positions - centroid

    # 結合の閾値（Å単位）
    bond_factor = float(data_toml["view"]["bond_factor"]) # 原子間距離がこの値以下なら結合とみなす
    bond_size_ratio = float(data_toml["view"]["bond_size_ratio"])
    atom_size_ratio = float(data_toml["view"]["atom_size_ratio"]) 
    
    # **結合リストを作る**
    distances = cdist(positions, positions)
    bonds = [(i, j) for i in range(len(positions)) for j in range(i+1, len(positions)) if distances[i, j] < bond_factor]

    data_atoms = Data_atoms()
    if not (data_toml["color_list"]["color_atoms"]):
        data_atoms.color[:] = ["#778899"] * len(data_atoms.color)
    
    # 色の設定
    #unique_numbers = list(set(numbers))
    #colors = {el: plt.cm.jet(i / len(unique_numbers)) for i, el in enumerate(unique_numbers)}
    
    # 球を描く関数
    def plot_sphere(ax, center, radius, color):
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 20)
        x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
        ax.plot_surface(x, y, z, color=color, alpha=1.0, edgecolor="k", linewidth=0.0, antialiased=True)

        vertices = np.array([x.ravel(), y.ravel(), z.ravel()]).T
        faces = []
        rows, cols = x.shape
        for i in range(rows - 1):
            for j in range(cols - 1):
                q1 = i * cols + j
                q2 = q1 + 1
                q3 = (i + 1) * cols + j
                q4 = q3 + 1
                faces.append([q1, q2, q3])  # 三角形1
                faces.append([q2, q4, q3])  # 三角形2

        return vertices, faces


    # **円柱（bond）をプロットする関数**
    def plot_cylinder(ax, p1, p2, radius=0.2, n=20, r_atom1=0.3, r_atom2=0.3):
        """ p1, p2 を結ぶ円柱を描画（原子半径を考慮） """
        v = np.array(p2) - np.array(p1)  # ベクトル p1 → p2
        length = np.linalg.norm(v)  # 円柱の高さ

        # **bond の端点を原子表面にシフト**
        delta_r1 = np.sqrt(float(r_atom1)*float(r_atom1) - float(radius)*float(radius))
        delta_r2 = np.sqrt(float(r_atom2)*float(r_atom2) - float(radius)*float(radius))
        p1_new = p1 + (delta_r1 / length) * v  # p1 から r_atom だけ v 方向へシフト
        p2_new = p2 - (delta_r2 / length) * v  # p2 から r_atom だけ v 方向へシフト
        v_new = p2_new - p1_new
        length_new = np.linalg.norm(v_new)  # 調整後の bond の長さ
    
        # 円柱の基本メッシュ (z 軸方向に伸びる円柱)
        theta = np.linspace(0, 2 * np.pi, n)
        z = np.linspace(0, length_new, 2)  # 円柱の高さ
        theta, z = np.meshgrid(theta, z)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
    
        # 方向を合わせる（回転）
        def rotation_matrix(v1, v2):
            """ v1 から v2 に回転する行列を求める """
            v1, v2 = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
            cross = np.cross(v1, v2)
            dot = np.dot(v1, v2)
            if np.allclose(cross, 0):  # 既に同じ方向なら回転不要
                return np.eye(3)
            skew = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
            return np.eye(3) + skew + (skew @ skew) * (1 / (1 + dot))
    
        # z軸基準の円柱を p1 → p2 に回転
        R = rotation_matrix(np.array([0, 0, 1]), v)
        xyz = np.dot(R, np.array([x.ravel(), y.ravel(), z.ravel()]))  # 座標変換
        x, y, z = xyz[0].reshape(x.shape) + p1_new[0], xyz[1].reshape(y.shape) + p1_new[1], xyz[2].reshape(z.shape) + p1_new[2]
    
        # プロット
        ax.plot_surface(x, y, z, color="white", alpha=1.0, linewidth=0.0, antialiased=True)

        vertices = np.array([x.ravel(), y.ravel(), z.ravel()]).T
        faces = []
        rows, cols = x.shape
        for i in range(rows - 1):
            for j in range(cols - 1):
                q1 = i * cols + j
                q2 = q1 + 1
                q3 = (i + 1) * cols + j
                q4 = q3 + 1
                faces.append([q1, q2, q3])  # 三角形1
                faces.append([q2, q4, q3])  # 三角形2

        return vertices, faces


    # **すべての結合を円柱で描画**
    for i, j in bonds:
        plot_cylinder(ax, positions[i], positions[j], radius=0.1*bond_size_ratio,
                                    r_atom1=data_atoms.size[numbers[i]]*atom_size_ratio,
                                    r_atom2=data_atoms.size[numbers[j]]*atom_size_ratio)

    # 各原子を球でプロット
    for pos, num in zip(positions, numbers):
        plot_sphere(ax, pos, data_atoms.size[num]*atom_size_ratio, data_atoms.color[num])

    ## **すべての結合を円柱で描画**
    #list_vertices_faces_bond = []
    #for i, j in bonds:
    #    list_vertices_faces_bond.append(plot_cylinder(ax, positions[i], positions[j], radius=0.1, r_atom=atom_radius))

    #mesh_list_bond = []
    #for tmp in list_vertices_faces_bond:
    #    mesh_list_bond.append(trimesh.Trimesh(vetices=tmp[0], faces=tmp[1]))


    ## 各原子を球でプロット
    #list_vertices_faces_atom= []
    #for pos, elem in zip(positions, numbers):
    #    list_vertices_faces_atom.append(plot_sphere(ax, pos, atom_radius, colors[elem]))

    #mesh_list_atom = []
    #for tmp in list_vertices_faces_atom:
    #    mesh_list_atom.append(trimesh.Trimesh(vetices=tmp[0], faces=tmp[1]))

    return


def create_supercell(verts, faces, supercell_size=(1, 1, 1)):
    """
    Marching CubesのvertsとfacesからSupercellを作成する関数
    
    Parameters:
    verts : numpy.ndarray
        (N, 3) の頂点座標リスト
    faces : numpy.ndarray
        (M, 3) の三角形のインデックス
    supercell_size : tuple
        (sx, sy, sz) の繰り返し数
    
    Returns:
    new_verts, new_faces : numpy.ndarray, numpy.ndarray
        Supercellの頂点と三角形インデックス
    """
    sx, sy, sz = supercell_size
    new_verts = []
    new_faces = []

    # Supercellの各シフトベクトルを計算
    shift_vectors = np.array([
        (i, j, k) for i in range(sx) for j in range(sy) for k in range(sz)
    ])

    num_original_verts = len(verts)
    
    for shift in shift_vectors:
        shift_vec = shift * np.array([1.0, 1.0, 1.0])  # 必要なら格子定数をかける
        new_verts.append(verts + shift_vec)

    new_verts = np.vstack(new_verts)

    # 各コピーごとにfacesを追加（インデックスはオフセットする）
    for i, shift in enumerate(shift_vectors):
        offset = i * num_original_verts
        new_faces.append(faces + offset)

    new_faces = np.vstack(new_faces)

    return new_verts, new_faces

def plot_isosurface(ax, data_toml, data_cube1, data_cube2):
    """ cube1の等値面を描き、その上にcube2のデータをカラーマップ＆等高線として描く """

    #print(repr(data_cube2))

    origin1       = data_cube1.origin
    grid_vectors1 = data_cube1.mesh_vec
    data1         = np.array(data_cube1.val).reshape(data_cube1.num_mesh)
    origin2       = data_cube2.origin
    grid_vectors2 = data_cube2.mesh_vec
    data2         = np.array(data_cube2.val).reshape(data_cube2.num_mesh)

    #origin1, grid_vectors1, data1 = read_cube(cube1)
    #origin2, grid_vectors2, data2 = read_cube(cube2)

    lattice1 = np.zeros((3,3))
    lattice1[0] = np.array(grid_vectors1[0]) * data_cube1.num_mesh[0]
    lattice1[1] = np.array(grid_vectors1[1]) * data_cube1.num_mesh[1]
    lattice1[2] = np.array(grid_vectors1[2]) * data_cube1.num_mesh[2]

    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
        # **重心を求める**
        centroid = np.mean(data_cube1.positions_ang, axis=0)

    # 3Dグリッド座標の作成
    #x = np.linspace(origin1[0], origin1[0] + grid_vectors1[0][0] * (data1.shape[0] - 1), data1.shape[0])
    #y = np.linspace(origin1[1], origin1[1] + grid_vectors1[1][1] * (data1.shape[1] - 1), data1.shape[1])
    #z = np.linspace(origin1[2], origin1[2] + grid_vectors1[2][2] * (data1.shape[2] - 1), data1.shape[2])
    #x2 = np.linspace(origin2[0], origin2[0] + grid_vectors2[0][0] * (data2.shape[0] - 1), data2.shape[0])
    #y2 = np.linspace(origin2[1], origin2[1] + grid_vectors2[1][1] * (data2.shape[1] - 1), data2.shape[1])
    #z2 = np.linspace(origin2[2], origin2[2] + grid_vectors2[2][2] * (data2.shape[2] - 1), data2.shape[2])

    from scipy.ndimage import gaussian_filter

    # ** 周期境界条件のために 1 周期分拡張 **
    wpad = [0,0,0]
    wpad[0] = int(data_toml["periodicity"][0])
    wpad[1] = int(data_toml["periodicity"][1])
    wpad[2] = int(data_toml["periodicity"][2])
    #data1_extended = np.pad(data1, pad_width=((1, 1), (1, 1), (1, 1)), mode='wrap')   
    data1_extended = np.pad(data1, pad_width=((wpad[0], wpad[0]), (wpad[1], wpad[1]), (wpad[2], wpad[2])), mode='wrap')   
    data2_extended = np.pad(data2, pad_width=((0, wpad[0]), (0, wpad[1]), (0, wpad[2])), mode='wrap')   

    x2_cubes = np.linspace(0, data_cube2.num_mesh[0], data_cube2.num_mesh[0]+wpad[0])
    y2_cubes = np.linspace(0, data_cube2.num_mesh[1], data_cube2.num_mesh[1]+wpad[1])
    z2_cubes = np.linspace(0, data_cube2.num_mesh[2], data_cube2.num_mesh[2]+wpad[2])

    if (data_toml["color_list"]["smooth"]):
        # データをスムージング
        smooth_data = gaussian_filter(data1_extended, sigma=0.8)  # sigmaを調整（大きいほど滑らか）

        # スムージング後のデータで等値面を抽出
        verts, faces, normals, _ = marching_cubes(smooth_data, level=data_toml["isovalue"], step_size=data_toml["view"]["roughness"])    
    else:
        # 等値面を抽出
        verts, faces, normals, _ = marching_cubes(data1_extended, level=data_toml["isovalue"], step_size=data_toml["view"]["roughness"])    

        #verts, faces2, _, _ = marching_cubes(data1, level=isovalue, step_size=3)
    
    # ** 元の範囲にある頂点のみをフィルタリング **
    valid_mask = np.logical_and.reduce([
        (verts[:, 0] >= wpad[0]) & (verts[:, 0] <= data1_extended.shape[0]-wpad[0]),
        (verts[:, 1] >= wpad[1]) & (verts[:, 1] <= data1_extended.shape[1]-wpad[1]),
        (verts[:, 2] >= wpad[2]) & (verts[:, 2] <= data1_extended.shape[2]-wpad[2])
    ])

    #valid_mask = np.logical_and.reduce([
    #    (verts[:, 0] >= 1) & (verts[:, 0] <= smooth_data.shape[0]-1),
    #    (verts[:, 1] >= 1) & (verts[:, 1] <= smooth_data.shape[1]-1),
    #    (verts[:, 2] >= 1) & (verts[:, 2] <= smooth_data.shape[2]-10)
    #])
      
    # mask後にフィルタリングしてずれたインデックス分をずらす。
    verts = verts[valid_mask]
    verts[:] = verts[:] - np.array([wpad[0],wpad[1],wpad[2]])

    # facesのフィルタリング
    face_mask = valid_mask[faces].all(axis=1)
    faces = faces[face_mask]
    
    # normalsのフィルタリング
    normals = normals[valid_mask]

    #normals = normals[valid_mask[faces].all(axis=1)]
    #faces = faces[valid_mask[faces].all(axis=1)]
    
    # valid_mask を使って、有効な頂点インデックスを再マッピングする
    valid_indices = np.where(valid_mask)[0]  # 有効な頂点のインデックスを取得
    
    # 高速リマッピング
    index_map = np.zeros(valid_mask.shape[0], dtype=int)
    index_map[valid_indices] = np.arange(valid_indices.shape[0])
    faces2 = index_map[faces]

    # 格子ベクトルで座標変換（ボクセル座標 → 実座標）
    verts_transformed = verts @ np.array(grid_vectors1)  
    verts_transformed[:] = verts_transformed[:] + origin1
    #verts_transformed = verts @ np.array(grid_vectors1).T  
    #x_surf_transformed = verts_transformed[:, 0].flatten()
    #y_surf_transformed = verts_transformed[:, 1].flatten()
    #z_surf_transformed = verts_transformed[:, 2].flatten()

    #verts_contour_transformed = verts_contour @ np.array(grid_vectors1)  
    #verts_contour_transformed[:] = verts_contour_transformed[:] + origin1
    #rx_s, ry_s, rz_s = verts_contour_transformed.T

    # 等値面の重心座標を取得
    x_surf, y_surf, z_surf = verts[faces2].mean(axis=1).T

    tmp_verts_transformed = copy.deepcopy(verts_transformed)
    tmp_verts_transformed[:] = tmp_verts_transformed[:] - origin2
    verts2 = tmp_verts_transformed @ np.linalg.inv(grid_vectors2)
    x_surf2, y_surf2, z_surf2 = verts2[faces2].mean(axis=1).T
    #print(verts)
    #print(x_surf.shape)
    #print(faces.shape)


    # cube2のデータを補間してカラーマッピング
    interp_func = RegularGridInterpolator((x2_cubes, y2_cubes, z2_cubes), data2_extended, bounds_error=False, fill_value=0)
    color_values = interp_func(np.c_[x_surf2, y_surf2, z_surf2])

    # しきい値 (s に近い点を抽出)
    #contour_level = -0.19
    #threshold = np.abs(contour_level/100)
    #selected_indices = np.where(np.abs(color_values - contour_level) < threshold)[0]
    #x_s = x_surf[selected_indices]
    #y_s = y_surf[selected_indices]
    #z_s = z_surf[selected_indices]
    #verts_contour = np.array([x_s, y_s, z_s]).T



    # ⭐ 数値データをRGBAカラーに変換する ⭐
    norm = colors.Normalize(vmin=data_toml["min_value"], vmax=data_toml["max_value"])
    #cmap = cm.get_cmap('RdYlBu')  # カラーマップを選択'coolwarm'
    cmap = colors.LinearSegmentedColormap.from_list(data_toml["color_list"]["name"], data_toml["color_list"]["colors"]) 
    face_colors_rgba = cmap(norm(color_values))  # RGBAに変換！

    #print(face_colors_rgba.shape)
    #print(faces.shape)
    #print(verts.shape)
    #print(verts[faces].shape)

#    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
#        # **重心を求める**
#        centroid = np.mean(verts_transformed, axis=0)
        
        # **全頂点を重心の逆方向に移動**
        verts_transformed = verts_transformed - centroid


    poly3d = Poly3DCollection(verts_transformed[faces2], facecolors=face_colors_rgba, linewidths=0.01, edgecolors=face_colors_rgba, alpha=data_toml["color_list"]["alpha"], shade=data_toml["color_list"]["shade"], zorder=10)
    ax.add_collection3d(poly3d)

    if (data_toml["show_lattice"]):
        draw_lattice_box(lattice1, origin1, centroid, ax)

    #rendering_by_blender(data_toml,verts_transformed, faces2)

    return


def plot_figure(data_toml):
    # 3Dプロット

    setting = Plot_style_1d()
    setting.plt_style()

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])  # X:Y:Z のスケールを等しくする

    # 軸ラベル
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ############  データのプロット  ##############
    data_cube1 = Data_cube()                                                                                                                                         
    data_cube1.read(data_toml["cube1"])
    data_cube2 = Data_cube()                                                                                                                                         
    data_cube2.read(data_toml["cube2"])

    if(data_toml["show_isosurface"]):
        plot_isosurface(ax, data_toml, data_cube1, data_cube2)
    if(data_toml["show_structure"]):
        plot_structure(ax, data_toml, data_cube1)

    ######################################


    ############  set plot range  ##############
    # 各軸の範囲を取得
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()
    
    # ★ 最大の範囲を求めて、すべての軸に適用
    #ratio_zoom = 1.5
    max_range = (max(x_max - x_min, y_max - y_min, z_max - z_min) / 2)/data_toml["view"]["ratio_zoom"]
    
    # ★ 各軸の中心を求める
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2
    
    # ★ 各軸の範囲を統一する
    ax.set_xlim([x_mid - max_range, x_mid + max_range])
    ax.set_ylim([y_mid - max_range, y_mid + max_range])
    ax.set_zlim([z_mid - max_range, z_mid + max_range])

    # 軸を非表示
    if not (data_toml["show_axis"]):
        ax.set_axis_off()


    if (data_toml["show_isosurface"]):
        if (data_toml["show_color_bar"]):
            ############  set color bar  ##############
            # カラーバーを追加
            norm = colors.Normalize(vmin=data_toml["min_value"], vmax=data_toml["max_value"])
            #cmap = cm.get_cmap('RdYlBu')  # カラーマップを選択'coolwarm'
            cmap = colors.LinearSegmentedColormap.from_list(data_toml["color_list"]["name"], data_toml["color_list"]["colors"]) 
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
            cbar.outline.set_linewidth(0.2)
            cbar.set_ticks([data_toml["min_value"], data_toml["max_value"]])
            #cbar.locator = MaxNLocator(nbins=2)  # 目盛りを5個に設定
            #cbar.update_ticks()  # 目盛りを更新
            # 目盛り線（tick marks）を非表示にする
            cbar.ax.tick_params(axis='both', which='both', length=0)


    ############  set view  ##############
    #ax.view_init(elev=00, azim=-90) #methane
    #ax.view_init(elev=-90, azim=120)
    ax.view_init(elev=data_toml["view"]["elev"], azim=data_toml["view"]["azim"])
    # 透視投影 (perspective) をオフにする（平行投影にする）
    ax.set_proj_type('ortho')


    if (data_toml["show_plot"]):
        #plt.show()
        plt.draw()
        plt.pause(1.0)

        # **vispy で OpenGL 描画**
        #from vispy import scene
        #from vispy.visuals import MeshVisual

        #canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', size=(600, 600))
        #view = canvas.central_widget.add_view()
        #mesh = MeshVisual(verts_transformed, faces2, color="red")
        ##mesh = MeshVisual(verts_transformed[faces2], faces2, face_colors=face_colors_rgba)

        ## `unfreeze()` を呼んで属性の追加を許可
        #mesh.unfreeze()

        #view.add(mesh)
        #view.camera = scene.cameras.TurntableCamera()  # マウスで回転
        #canvas.show()


    if (data_toml["save_pdf"]):
        plt.savefig(data_toml["output_name"] + ".pdf", bbox_inches='tight')
        #plt.savefig("chempot.pdf", bbox_inches='tight', pad_inches=0)

    if (data_toml["save_obj"]):
        write_obj_and_mtl(verts_transformed, faces2, face_colors_rgba, \
                          obj_filename=data_toml["output_name"] + ".obj", \
                          mtl_filename=data_toml["output_name"] + ".mtl")
        #export_gltf(verts_transformed.astype(np.float32), face_colors_rgba.astype(np.float32), faces2.astype(np.int32), normals.astype(np.float32))


def write_obj_and_mtl(vertices, faces, colors, obj_filename="model.obj", mtl_filename="model.mtl"):
    # MTLファイルを作成
    with open(mtl_filename, "w") as mtl_file:
        mtl_file.write("# MTL file\n")

        # マテリアルを作成
        for i, color in enumerate(colors):
            material_name = f"Material_{i+1}"
            mtl_file.write(f"newmtl {material_name}\n")
            mtl_file.write(f"Ka {color[0]} {color[1]} {color[2]}\n")  # 環境光の色
            mtl_file.write(f"Kd {color[0]} {color[1]} {color[2]}\n")  # 拡散光の色
            mtl_file.write(f"Ks 1.0 1.0 1.0\n")  # 鏡面反射の色（白）
            mtl_file.write("\n")

    # OBJファイルを作成
    with open(obj_filename, "w") as obj_file:
        obj_file.write("# OBJ file\n")
        obj_file.write(f"mtllib {mtl_filename}\n")  # MTLファイルを参照

        # 頂点の書き込み
        for v in vertices:
            obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")

        # 面の書き込み (OBJは1-based index)
        for i, face in enumerate(faces):
            material_idx = i  # 各面に色を割り当てる
            obj_file.write(f"usemtl Material_{material_idx+1}\n")  # マテリアルを設定
            obj_file.write(f"f {' '.join(str(i+1) for i in face)}\n")


def draw_lattice_box(lattice, origin, centroid, ax):
    
    # ボックスの8つの頂点を計算
    v1 = lattice[0]
    v2 = lattice[1] 
    v3 = lattice[2]
    
    box_vertices = np.array([
        np.zeros(3),        
        v1, v2, v3, 
        v1 + v2, v1 + v3, v2 + v3,
        v1 + v2 + v3
    ])

    box_vertices = box_vertices + np.array(origin) - centroid
    
    # ボックスの12本のエッジを定義
    edges = [
        [0, 1], [0, 2], [0, 3],  # 原点からのエッジ
        [1, 4], [1, 5],          # v1 からのエッジ
        [2, 4], [2, 6],          # v2 からのエッジ
        [3, 5], [3, 6],          # v3 からのエッジ
        [4, 7], [5, 7], [6, 7]   # 最後の3つのエッジ
    ]

    #for edge in edges:
    #    p1, p2 = box_vertices[edge]
    #    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 'k', lw=0.5, alpha=1.0, zorder=10)
    
    # ** 
    # 視点の方向ベクトルを計算
    elev=data_toml["view"]["elev"]
    azim=data_toml["view"]["azim"]
    view_vector = np.array([np.cos(np.radians(azim)) * np.cos(np.radians(elev)),
                        np.sin(np.radians(azim)) * np.cos(np.radians(elev)),
                        np.sin(np.radians(elev))])

    # 各線の中心点の Z値(深度)を計算してソート
    z_values = np.array([np.dot(box_vertice, view_vector) for box_vertice in box_vertices])
    sorted_indices = np.argsort(z_values)  # 手前から順に並べる


    # ボックスのエッジを描画
    lines_front = []
    lines_behind = []
    for edge in edges:
        p1, p2 = box_vertices[edge]
        if((edge[0]==sorted_indices[0]) or (edge[1]==sorted_indices[0])):
            lines_front.append([p1,p2])
        else:
            lines_behind.append([p1,p2])

    
    # **Line3DCollection の作成**
    line_collection_front = Line3DCollection(lines_front, colors="black", linewidths=0.5, zorder=1)
    line_collection_behind = Line3DCollection(lines_behind, colors="black", linewidths=0.5, zorder=100)
    
    # **プロットに追加**
    ax.add_collection3d(line_collection_front)
    ax.add_collection3d(line_collection_behind)

    ## ボックスの面を半透明で描画（オプション）
    #faces = [[box_vertices[i] for i in face] for face in [
    #    [0, 1, 4, 2], [0, 1, 5, 3], [0, 2, 6, 3],
    #    [1, 4, 7, 5], [2, 4, 7, 6], [3, 5, 7, 6]
    #]]
    #ax.add_collection3d(Poly3DCollection(faces, alpha=0.3, facecolor='cyan'))
    

def map_colors(p3dc, func, cmap='viridis'):
    """
Color a tri-mesh according to a function evaluated in each barycentre.

    p3dc: a Poly3DCollection, as returned e.g. by ax.plot_trisurf
    func: a single-valued function of 3 arrays: x, y, z
    cmap: a colormap NAME, as a string

    Returns a ScalarMappable that can be used to instantiate a colorbar.
    """
    
    from matplotlib.cm import ScalarMappable, get_cmap
    from matplotlib.colors import Normalize
    from numpy import array

    # reconstruct the triangles from internal data
    x, y, z, _ = p3dc._vec
    slices = p3dc._segslices
    triangles = array([array((x[s],y[s],z[s])).T for s in slices])

    # compute the barycentres for each triangle
    xb, yb, zb = triangles.mean(axis=1).T
    
    # compute the function in the barycentres
    values = func(xb, yb, zb)

    # usual stuff
    norm = Normalize()
    colors = get_cmap(cmap)(norm(values))

    # set the face colors of the Poly3DCollection
    p3dc.set_fc(colors)

    # if the caller wants a colorbar, they need this
    return ScalarMappable(cmap=cmap, norm=norm)

def export_gltf(vertices, colors, faces, normals):

    import struct
    import json

    print(vertices.dtype)
    print(colors.dtype)
    print(faces.dtype)
    print(normals.dtype)

    # バイナリバッファ作成
    bin_data = vertices.tobytes() + colors.tobytes() + faces.tobytes()
    #bin_data = vertices.tobytes() + normals.tobytes() + colors.tobytes() + faces.tobytes()

    gltf_json = {
        "asset": {"version": "2.0"},
        "scenes": [{"nodes": [0]}],
        "nodes": [{"mesh": 0}],
        "meshes": [{
            "primitives": [{
                "attributes": {
                    "POSITION": 0,
                    "COLOR_0": 1
                },
                "indices": 2
            }]
        }],
        "accessors": [
            {
                "bufferView": 0, "componentType": 5126, "count": len(vertices),
                "type": "VEC3", "max": vertices.max(axis=0).tolist(), "min": vertices.min(axis=0).tolist()
            },
            {
                "bufferView": 1, "componentType": 5126, "count": len(colors),
                "type": "VEC4"
            },
            {
                "bufferView": 2, "componentType": 5125, "count": len(faces) * 3,
                "type": "SCALAR"
            }
        ],
        "bufferViews": [
            {"buffer": 0, "byteOffset": 0, "byteLength": vertices.nbytes},
            {"buffer": 0, "byteOffset": vertices.nbytes, "byteLength": colors.nbytes},
            {"buffer": 0, "byteOffset": vertices.nbytes + colors.nbytes, "byteLength": faces.nbytes}
        ],
        "buffers": [{"byteLength": vertices.nbytes + colors.nbytes + faces.nbytes}]
    }

    # JSON データ（GLBヘッダー）
    #gltf_json = {
    #    "asset": {"version": "2.0"},
    #    "scenes": [{"nodes": [0]}],
    #    "nodes": [{"mesh": 0}],
    #    "meshes": [{
    #        "primitives": [{
    #            "attributes": {
    #                "POSITION": 0,
    #                "NORMAL": 1,  # 法線のインデックスを修正
    #                "COLOR_0": 2
    #            },
    #            "indices": 3  # 面のインデックスも修正
    #        }]
    #    }],
    #    "accessors": [
    #        {
    #            "bufferView": 0, "componentType": 5126, "count": len(vertices),
    #            "type": "VEC3", "max": vertices.max(axis=0).tolist(), "min": vertices.min(axis=0).tolist()
    #        },
    #        {
    #            "bufferView": 1, "componentType": 5126, "count": len(normals),
    #            "type": "VEC3"
    #        },
    #        {
    #            "bufferView": 2, "componentType": 5126, "count": len(colors),
    #            "type": "VEC4"
    #        },
    #        {
    #            "bufferView": 3, "componentType": 5125, "count": len(faces) * 3,
    #            "type": "SCALAR"
    #        }
    #    ],
    #    "bufferViews": [
    #        {"buffer": 0, "byteOffset": 0, "byteLength": vertices.nbytes},  # 頂点
    #        {"buffer": 0, "byteOffset": vertices.nbytes, "byteLength": normals.nbytes},  # 法線
    #        {"buffer": 0, "byteOffset": vertices.nbytes + normals.nbytes, "byteLength": colors.nbytes},  # 色
    #        {"buffer": 0, "byteOffset": vertices.nbytes + normals.nbytes + colors.nbytes, "byteLength": faces.nbytes}  # 面
    #    ],
    #    "buffers": [{"byteLength": len(bin_data)}]
    #}

    # JSON 保存
    with open("model.gltf", "w") as f:
        json.dump(gltf_json, f)

    print(len(bin_data))  # 実際のバイナリデータのサイズ
    print(gltf_json["buffers"][0]["byteLength"])  # JSON に書いたサイズ
    
    # JSON をバイナリに変換
    json_data = json.dumps(gltf_json).encode('utf-8')
    json_length = len(json_data)
    bin_length = len(bin_data)
    
    # GLB ヘッダー
    glb_header = struct.pack("<I", 0x46546C67)  # "glTF" マジックナンバー
    glb_header += struct.pack("<I", 2)          # GLTF バージョン
    glb_header += struct.pack("<I", 12 + 8 + json_length + 8 + bin_length)  # ファイルサイズ
    
    # JSON チャンク
    json_chunk_header = struct.pack("<I", json_length)
    json_chunk_header += struct.pack("<I", 0x4E4F534A)  # "JSON"
    json_chunk_data = json_data
    
    # バイナリ チャンク
    bin_chunk_header = struct.pack("<I", bin_length)
    bin_chunk_header += struct.pack("<I", 0x004E4942)  # "BIN"
    bin_chunk_data = bin_data
    
    # GLB 書き込み
    with open("model.glb", "wb") as f:
        f.write(glb_header)
        f.write(json_chunk_header)
        f.write(json_chunk_data)
        f.write(bin_chunk_header)
        f.write(bin_chunk_data)

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

def rendering_by_blender(data_toml,vertices, faces):

    import bpy

    #vertices_N3 = np.column_stack([v.ravel() for v in vertices])

    # NumPy → Blender 用にタプルのリストに変換
    #vertices_N3 = [tuple(v) for v in vertices_N3.tolist()]

    #print(vertices.dtype)
    #print(vertices.shape)  # 配列の形状を確認
    #print(vertices[:5])    # 最初の5つを表示

    # **既存のオブジェクトを削除**
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

    ## **OBJ ファイルをインポート**
    #bpy.ops.import_scene.obj(filepath=data_toml["output_name"] + ".obj")
    #obj = bpy.context.selected_objects[0]

    # **新しいメッシュを作成**
    mesh = bpy.data.meshes.new(name="MyMesh")
    obj = bpy.data.objects.new("MyObject", mesh)
    
    # **シーンに追加**
    bpy.context.collection.objects.link(obj)
    
    # **メッシュデータを設定**
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    
    # **スムーズシェーディング**
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.ops.object.shade_smooth()

    # **マテリアルを作成**
    mat = bpy.data.materials.new(name="MyMaterial")
    mat.diffuse_color = (1.0, 0.5, 0.2, 1.0)  # オレンジ色
    obj.data.materials.append(mat)
    
    # **ライトを追加**
    bpy.ops.object.light_add(type='SUN', location=(5, -5, 10))
    bpy.context.object.data.energy = 5  # 明るさ
    # 追加されたライトを取得
    sun_light = bpy.context.object
    # `SUN` ライトの影を無効にする
    sun_light.data.use_shadow = False


    r = 10     # カメラの距離

    # **角度をラジアンに変換**
    elev_rad = np.radians(data_toml["view"]["elev"])
    azim_rad = np.radians(data_toml["view"]["azim"])
    
    # **カメラ位置を計算**
    x = r * np.cos(azim_rad) * np.cos(elev_rad)
    y = r * np.sin(azim_rad) * np.cos(elev_rad)
    z = r * np.sin(elev_rad)
    
    # **Blender でカメラを設定**
    bpy.ops.object.camera_add()
    camera = bpy.context.object
    camera.location = (x, y, z)
    
    # **シーンのメインカメラとして登録**
    bpy.context.scene.camera = camera
    
    # **ターゲット用の "Empty" オブジェクトを作成（原点に設置）**
    bpy.ops.object.empty_add(location=(0, 0, 0))
    target = bpy.context.object
    
    # **カメラの向きをターゲットに向ける**
    camera_constraint = camera.constraints.new(type='TRACK_TO')
    camera_constraint.target = target
    camera_constraint.track_axis = 'TRACK_NEGATIVE_Z'
    camera_constraint.up_axis = 'UP_Y'

    # **レンダリング設定**
    bpy.context.scene.render.resolution_x = 1920  # 横幅
    bpy.context.scene.render.resolution_y = 1080  # 縦幅
    bpy.context.scene.render.resolution_percentage = 50  # スケール (100%でフル解像度)

    bpy.context.scene.render.engine = 'CYCLES'  # 高品質レンダリング
    bpy.context.scene.render.image_settings.file_format = 'PNG'
    output_dir = os.path.dirname(bpy.data.filepath) if bpy.data.filepath else os.getcwd()
    output_path = os.path.join(output_dir, "output.png")
    bpy.context.scene.render.filepath = output_path  # 保存パス
    
    # **レンダリング開始**
    bpy.ops.render.render(write_still=True)


def get_bibtex():
    bibtex_entries = {
        "trimesh": """@article{trimesh,
    author  = {Dawson-Haggerty, M. and others},
    title   = {trimesh: Python library for loading and using triangular meshes},
    year    = {2019},
    note    = {Available at \\url{https://github.com/mikedh/trimesh}}
}""",
        "scikit-image": """@article{scikit-image,
    author  = {Van der Walt, S. and Sch\\"{o}nberger, J. L. and Nunez-Iglesias, J. and Boulogne, F. and Warner, J. D. and others},
    title   = {scikit-image: image processing in Python},
    journal = {PeerJ},
    volume  = {2},
    pages   = {e453},
    year    = {2014},
    doi     = {10.7717/peerj.453}
}""",
        "scipy": """@article{scipy,
    author  = {Virtanen, P. and Gommers, R. and Oliphant, T. E. and Haberland, M. and Reddy, T. and others},
    title   = {{SciPy} 1.0: fundamental algorithms for scientific computing in Python},
    journal = {Nature Methods},
    volume  = {17},
    pages   = {261--272},
    year    = {2020},
    doi     = {10.1038/s41592-019-0686-2}
}"""
    }
    return "\n\n".join(bibtex_entries.values())

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
    def __init__(self, toml_file, default_config, data_toml):
        super().__init__()
        self.toml_file = toml_file
        self.default_config = default_config

    def on_modified(self, event):
        global check_update_flag
        if event.src_path.endswith(self.toml_file):
            if (check_update_flag) :
                user_config = read_toml(self.toml_file)
                data_toml = deep_update(self.default_config, user_config)
                print("toml file is updated")

                plot_queue.put(data_toml)
                check_update_flag = False
            else:
                check_update_flag = True


# 監視を開始する関数
def start_watching(toml_file, default_config, data_toml):

    event_handler = TomlChangeHandler(toml_file, default_config, data_toml)
    observer = Observer()
    observer.schedule(event_handler, ".", recursive=False)
    observer.start()
    
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

def plot_figure_interactive():
    while True:
        try:
            # キューからデータを取得
            data_toml = plot_queue.get(timeout=1)
            plt.close("all")
            plot_figure(data_toml)

        except queue.Empty:
            continue  # キューが空ならスキップ

if __name__ == '__main__':

    import sys

    parser = argparse.ArgumentParser(description="Generate BibTeX for cited Python modules.")
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


    if(data_toml["show_plot"]):
        # メインスレッドで実行するためのキュー
        plot_queue = queue.Queue()

        plot_queue.put(data_toml)

        # 別スレッドで監視を開始
        thread = threading.Thread(target=start_watching, args=(toml_file, default_config, data_toml), daemon=True)
        thread.start()
        
        # メインスレッドでプロット更新を監視するスレッドを開始
        plot_thread = threading.Thread(target=plot_figure_interactive(), daemon=True)
        plot_thread.start()


        # メインスレッドをブロック（終了しないように）
        while True:
            time.sleep(1)


    plot_figure(data_toml)

    #plot_isosurface("LPQ.kinetic_energy_density.cube", "LPQ.chemical_potential.cube", max_value=-0.20, min_value=-0.38, isovalue=0.00001, contour_levels=5)
    #plot_isosurface("LPQ.kinetic_energy_density.cube", "LPQ.chemical_potential.cube", max_value=-0.175, min_value=-0.2, isovalue=0.00001, contour_levels=5)
    #plot_isosurface("RESULT.tden.cube", "RESULT.vhart.cube", max_value=0.02, min_value=-0.02, isovalue=0.0005)

