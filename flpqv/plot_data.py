import os
import pprint
import copy
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


#
from parameters import BohrR
from parameters import Data_atoms


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

    #rendering_by_blender(data_toml,verts_transformed, faces2)

    return


def plot_lattice_box(ax, data_toml, data_cube):

    origin = data_cube.origin

    lattice = np.zeros((3,3))
    lattice[0] = np.array(data_cube.mesh_vec[0]) * data_cube.num_mesh[0]
    lattice[1] = np.array(data_cube.mesh_vec[1]) * data_cube.num_mesh[1]
    lattice[2] = np.array(data_cube.mesh_vec[2]) * data_cube.num_mesh[2]

    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
        # **重心を求める**
        centroid = np.mean(data_cube.positions_ang, axis=0)
    
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


