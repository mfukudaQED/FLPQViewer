import os
import pprint
import copy
import numpy as np

import trimesh

from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotly.graph_objects as go



#
from parameters import BohrR
from parameters import Data_atoms

def hex_to_rgb(hex_code):
    hex_code = hex_code.lstrip("#")  # "#" を削除
    return tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))
    # print(hex_to_rgb("#FF5733"))  # (255, 87, 51)

def get_text_color(rgb):
    """Determine whether to use black or white text based on background color."""
    rgb01 = np.array(rgb)/255
    luminance = 0.299 * rgb01[0] + 0.587 * rgb01[1] + 0.114 * rgb01[2]
    return "black" if luminance > 0.5 else "white"

def plot_structure(figp, ax, text_objects, data_toml, data_cube):
    """
    原子構造をプロットするコード
    """

    toml_s = data_toml["structure"]

    # 原子座標を取得
    positions = data_cube.positions_ang
    numbers = np.int64(data_cube.numbers)
    indices = data_cube.indices

    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
        # **重心を求める**
        centroid = np.mean(data_cube.positions_ang, axis=0)

    positions = positions - centroid + np.array(data_toml["view"]["translate"])

    # 結合の閾値（Å単位）
    bond_factor = float(toml_s["bond_factor"]) # 原子間距離がこの値以下なら結合とみなす
    bond_size_ratio = float(toml_s["bond_size_ratio"])
    atom_size_ratio = float(toml_s["atom_size_ratio"]) 
    
    # **結合リストを作る**
    distances = cdist(positions, positions)
    bonds = [(i, j) for i in range(len(positions)) for j in range(i+1, len(positions)) if distances[i, j] < bond_factor]

    data_atoms = Data_atoms()
    if not (toml_s["color_atoms"]):
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
        ax.plot_surface(x, y, z, color=color, alpha=1.0, edgecolor="k", linewidth=0.0, antialiased=True, zorder=25)

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
    def plot_cylinder(ax, norm, cmap, p1, p2, radius=0.2, n=20, r_atom1=0.3, r_atom2=0.3):
        """ p1, p2 を結ぶ円柱を描画（原子半径を考慮） """
        v = np.array(p2) - np.array(p1)  # ベクトル p1 → p2
        length = np.linalg.norm(v)  # 円柱の高さ
        if (toml_s["color_bonds"]):
            color = cmap(norm(length))  # RGBAに変換！
        else:
            color = "white"

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
        ax.plot_surface(x, y, z, color=color, alpha=1.0, linewidth=0.0, antialiased=True, zorder=25)

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

        return vertices, faces, color, length



    norm = colors.Normalize(vmin=toml_s["min_bond"], vmax=toml_s["max_bond"])
    cmap = colors.LinearSegmentedColormap.from_list(toml_s["color_name_bond"], toml_s["colors_bond"]) 
    if (toml_s["color_bar_bond"]):
        ############  set color bar  ##############
        # カラーバーを追加
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        #divider = make_axes_locatable(ax)
        #cax_bond = divider.append_axes("right", size="5%", pad=0.1)
        #cbar = plt.colorbar(sm, ax=cax_bond, shrink=0.5, aspect=20, label="bond")
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20, label="bond")
        cbar.outline.set_linewidth(0.2)
        cbar.set_ticks([toml_s["min_bond"], toml_s["max_bond"]])
        cbar.locator = MaxNLocator(nbins=2)  # 目盛りを5個に設定
        #cbar.update_ticks()  # 目盛りを更新
        # 目盛り線（tick marks）を非表示にする
        cbar.ax.tick_params(axis='both', which='both', length=0)

    # **すべての結合を円柱で描画**
    list_vertices_faces_bond= []
    for i, j in bonds:
        vertices, faces, color, length = plot_cylinder(ax, norm, cmap, positions[i], positions[j], radius=0.1*bond_size_ratio,
                                                r_atom1=data_atoms.size[numbers[i]]*atom_size_ratio,
                                                r_atom2=data_atoms.size[numbers[j]]*atom_size_ratio)
        list_vertices_faces_bond.append((vertices, faces, color, length))


    #proj = ax.get_proj()  # Get the projection matrix
    #scale = max(np.abs(proj[0, 0]), np.abs(proj[1, 1]), np.abs(proj[2, 2]))

    # 各原子を球でプロット
    list_vertices_faces_atom= []
    list_face_colors_atom = []
    for index, pos, num in zip(indices, positions, numbers):
        list_vertices_faces_atom.append(plot_sphere(ax, pos, data_atoms.size[num]*atom_size_ratio, data_atoms.color[num]))
        list_face_colors_atom.append(data_atoms.color[num]) 
        text_color = get_text_color(hex_to_rgb(data_atoms.color[num]))
        #print(text_color)
        text_fontsize = 24*data_atoms.size[num]*atom_size_ratio*data_toml["view"]["ratio_zoom"]
        #text_fontsize_scaled =text_fontsize / scale
        if (toml_s["number"]):
            text = ax.text(pos[0], pos[1], pos[2], index, color=text_color, fontsize=text_fontsize, ha='center', va='center', zorder=75)
            text_objects.append((text, text_color, text_fontsize))
        elif (toml_s["element"]):
            text = ax.text(pos[0], pos[1], pos[2], data_atoms.name[num], color=text_color, fontsize=text_fontsize, ha='center', va='center', zorder=75)
            text_objects.append((text, text_color, text_fontsize))

    #list_vertices_atom = np.vstack(list_vertices_faces_atom[0])
    #list_faces_atom = np.vstack(list_vertices_faces_atom[1])

    #layout = go.Layout(
    #    scene=dict(
    #        xaxis=dict(title='X', range=[-5, 5]),  # Set range for X axis
    #        yaxis=dict(title='Y', range=[-5, 5]),  # Set range for Y axis
    #        zaxis=dict(title='Z', range=[-5, 5]),  # Set range for Z axis
    #        camera=dict(
    #            projection=dict(
    #                type='orthographic'  # Use orthographic projection
    #            ),
    #            eye=dict(x=1.5, y=1.5, z=1.5)  # Adjust camera position
    #        )
    #    ),
    #)

    #figp = go.Figure(layout=layout)
    for index, tmp in enumerate(list_vertices_faces_atom):
        i, j, k = zip(*tmp[1])

        figp.add_trace(go.Mesh3d(
            x=tmp[0][:, 0],  # x座標
            y=tmp[0][:, 1],  # y座標
            z=tmp[0][:, 2],  # z座標
            i=i,     # 面の1番目の頂点インデックス
            j=j,     # 面の2番目の頂点インデックス
            k=k,     # 面の3番目の頂点インデックス
            color = list_face_colors_atom[index],
            opacity=1.0,       # 透明度
            flatshading=False,  # フラットシェーディングを無効にする
            lighting=dict(
                ambient=0.7,   # 環境光（全体の明るさ）
                diffuse=0.9,   # 拡散光（物体の明るさ）
                specular=1.0,  # 鏡面反射（光の反射）
                roughness=0.2,  # 表面の粗さ
                fresnel=0.2     # フレネル効果
            ),
            #lightposition=dict(
            #    x=100,  # 光の位置（X軸）
            #    y=100,  # 光の位置（Y軸）
            #    z=100   # 光の位置（Z軸）
            #)
        ))

    for tmp in list_vertices_faces_bond:
        i, j, k = zip(*tmp[1])
        if (toml_s["color_bonds"]):
            color = f"rgb({255*tmp[2][0]}, {255*tmp[2][1]}, {255*tmp[2][2]})"
        else:
            color = "white"


        figp.add_trace(go.Mesh3d(
            x=tmp[0][:, 0],  # x座標
            y=tmp[0][:, 1],  # y座標
            z=tmp[0][:, 2],  # z座標
            i=i,     # 面の1番目の頂点インデックス
            j=j,     # 面の2番目の頂点インデックス
            k=k,     # 面の3番目の頂点インデックス
            color = color,
            opacity=1.0,       # 透明度
            flatshading=False,  # フラットシェーディングを無効にする
            lighting=dict(
                ambient=0.7,   # 環境光（全体の明るさ）
                diffuse=0.9,   # 拡散光（物体の明るさ）
                specular=1.0,  # 鏡面反射（光の反射）
                roughness=0.2,  # 表面の粗さ
                fresnel=0.2     # フレネル効果
            ),
            #lightposition=dict(
            #    x=100,  # 光の位置（X軸）
            #    y=100,  # 光の位置（Y軸）
            #    z=100   # 光の位置（Z軸）
            #)
        ))

    #figp.write_html("output.html")

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


def create_supercell(lattice, verts, faces=None, face_colors=None, supercell_size=(1, 1, 1)):
    """
    Create a supercell from the given vertices and faces.

    Parameters:
    lattice : numpy.ndarray
        (3, 3) lattice matrix defining the unit cell
    verts : numpy.ndarray
        (N, 3) array of vertex coordinates
    faces : numpy.ndarray
        (M, 3) array of triangle indices
    face_colors : numpy.ndarray
        (M, 3) array of colors corresponding to each face
    supercell_size : tuple
        (sx, sy, sz) defining the repetition count along each axis

    Returns:
    new_lattice : numpy.ndarray
        (3, 3) lattice matrix of the supercell
    new_verts : numpy.ndarray
        (N * sx * sy * sz, 3) array of vertex coordinates
    new_faces : numpy.ndarray
        (M * sx * sy * sz, 3) array of triangle indices
    new_face_colors : numpy.ndarray
        (M * sx * sy * sz, 3) array of colors for the supercell faces
    """

    sx, sy, sz = supercell_size
    if (sx==1 and sy==1 and sz==1):
        return lattice, verts, faces, face_colors

    new_verts = []
    new_faces = [] if faces is not None else None
    new_face_colors = [] if face_colors is not None else None

    # Supercellの各シフトベクトルを計算
    shift_vectors = np.array([
        (i, j, k) for i in range(sx) for j in range(sy) for k in range(sz)
    ])

    num_original_verts = len(verts)
    
    for shift in shift_vectors:
        shift_vec = np.dot(shift, lattice)
        new_verts.append(verts + shift_vec)

    new_verts = np.vstack(new_verts)

    # 各コピーごとにfacesを追加（インデックスはオフセットする）
    if faces is not None:
        for i, shift in enumerate(shift_vectors):
            offset = i * num_original_verts
            new_faces.append(faces + offset)
            if face_colors is not None:
                new_face_colors.append(face_colors)

    if faces is not None:
        new_faces = np.vstack(new_faces)
        if face_colors is not None:
            new_face_colors = np.vstack(new_face_colors)

    # Scale lattice for supercell
    new_lattice = lattice * np.array(supercell_size)[:, np.newaxis] 

    return new_lattice, new_verts, new_faces, new_face_colors

def plot_isosurface(figp, ax, data_toml, data_cube1, data_cube2, isovalue):
    """ cube1の等値面を描き、その上にcube2のデータをカラーマップ＆等高線として描く """

    toml_i = data_toml["isosurface"]

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
    wpad[0] = int(toml_i["periodicity"][0])
    wpad[1] = int(toml_i["periodicity"][1])
    wpad[2] = int(toml_i["periodicity"][2])
    #data1_extended = np.pad(data1, pad_width=((1, 1), (1, 1), (1, 1)), mode='wrap')   
    data1_extended = np.pad(data1, pad_width=((wpad[0], wpad[0]), (wpad[1], wpad[1]), (wpad[2], wpad[2])), mode='wrap')   
    data2_extended = np.pad(data2, pad_width=((0, wpad[0]), (0, wpad[1]), (0, wpad[2])), mode='wrap')   

    x2_cubes = np.linspace(0, data_cube2.num_mesh[0], data_cube2.num_mesh[0]+wpad[0])
    y2_cubes = np.linspace(0, data_cube2.num_mesh[1], data_cube2.num_mesh[1]+wpad[1])
    z2_cubes = np.linspace(0, data_cube2.num_mesh[2], data_cube2.num_mesh[2]+wpad[2])

    if (toml_i["smooth"]):
        # データをスムージング
        smooth_data = gaussian_filter(data1_extended, sigma=0.8)  # sigmaを調整（大きいほど滑らか）

        # スムージング後のデータで等値面を抽出
        #verts, faces, normals, _ = marching_cubes(smooth_data, level=toml_i["isovalue"], step_size=data_toml["view"]["roughness"])    
        verts, faces, normals, _ = marching_cubes(smooth_data, level=isovalue, step_size=data_toml["view"]["roughness"])    

        valid_mask = np.logical_and.reduce([
            (verts[:, 0] >= wpad[0]) & (verts[:, 0] <= smooth_data.shape[0]-wpad[0]),
            (verts[:, 1] >= wpad[1]) & (verts[:, 1] <= smooth_data.shape[1]-wpad[1]),
            #(verts[:, 2] >= wpad[2]) & (verts[:, 2] <= smooth_data.shape[2]-wpad[2])
            (verts[:, 2] >= wpad[2]) & (verts[:, 2] <= smooth_data.shape[2]-wpad[2]-23)
        ])

    else:
        # 等値面を抽出
        verts, faces, normals, _ = marching_cubes(data1_extended, level=isovalue, step_size=data_toml["view"]["roughness"])    

        #verts, faces2, _, _ = marching_cubes(data1, level=isovalue, step_size=3)
    
        # ** 元の範囲にある頂点のみをフィルタリング **
        valid_mask = np.logical_and.reduce([
            (verts[:, 0] >= wpad[0]) & (verts[:, 0] <= data1_extended.shape[0]-wpad[0]),
            (verts[:, 1] >= wpad[1]) & (verts[:, 1] <= data1_extended.shape[1]-wpad[1]),
            (verts[:, 2] >= wpad[2]) & (verts[:, 2] <= data1_extended.shape[2]-wpad[2])
        ])

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
    #color_values = interp_func(np.c_[x_surf2, y_surf2, z_surf2])


    # すべての頂点のz座標を集める
    all_z_values = [vertex[2] for face in faces2 for vertex in verts2[face]]
    
    # 最大値を計算
    #max_z = max(all_z_values)
    min_z = min(all_z_values)
    max_z = max(all_z_values)

    if (toml_i["plot_2d"]):
        # Height of plane
        #z_plane = min_z - 1.5
        #z_plane = min_z - toml_i["plot_2d_z"]
        z_plane = max_z + toml_i["plot_2d_z"]
        r_s = 4.0
        print(f"max_z: {max_z}")
        print(f"plot_2d_z: {toml_i["plot_2d_z"]}")
        print(f"z_plane: {z_plane}")
        #alpha = 0.3

    color_values = []
    for face in faces2:
        tmp = np.float64(0.0)
        for vertex in verts2[face]:
            if (toml_i["plot_2d"]):
                tmp_mu = interp_func(np.c_[vertex[0], vertex[1], vertex[2]])
                #alpha = np.sqrt(np.abs(tmp_mu))
                alpha = -tmp_mu*2.0*r_s
                #print(f"alpha: {alpha}")
                #print(f"z: {vertex[2]}")
                tmp = tmp - (1.0-alpha*(z_plane - vertex[2] + r_s)) *np.exp( -alpha*( z_plane - vertex[2] + r_s))
                #tmp = tmp + (1.0+tmp_mu) *np.exp( -alpha*( z_plane - vertex[2] ))
                ##tmp = tmp + tmp_mu *np.exp( alpha*( z_plane - vertex[2] ))
                ##tmp = tmp - np.sqrt(np.abs(tmp_mu)) *np.exp( -alpha*( z_plane - vertex[2] ))
                ##tmp = tmp - np.sqrt(np.abs(tmp_mu)) *np.exp( alpha*( z_plane - vertex[2] ))

            else:
                tmp = tmp + interp_func(np.c_[vertex[0], vertex[1], vertex[2]])

        tmp = tmp/len(verts2[face])
        color_values.append(tmp)

    # しきい値 (s に近い点を抽出)
    #contour_level = -0.19
    #threshold = np.abs(contour_level/100)
    #selected_indices = np.where(np.abs(color_values - contour_level) < threshold)[0]
    #x_s = x_surf[selected_indices]
    #y_s = y_surf[selected_indices]
    #z_s = z_surf[selected_indices]
    #verts_contour = np.array([x_s, y_s, z_s]).T



    # ⭐ 数値データをRGBAカラーに変換する ⭐
    norm = colors.Normalize(vmin=toml_i["min_value"], vmax=toml_i["max_value"])
    #cmap = cm.get_cmap('RdYlBu')  # カラーマップを選択'coolwarm'
    cmap = colors.LinearSegmentedColormap.from_list(toml_i["color_name"], toml_i["colors"]) 
    face_colors_rgba = cmap(norm(color_values))  # RGBAに変換！

    lattice1, verts_transformed, faces2, face_colors_rgba = create_supercell(lattice1, verts_transformed, faces2, face_colors_rgba, supercell_size=data_toml["view"]["supercell"])

    if (toml_i["color_bar"]):
        ############  set color bar  ##############
        # カラーバーを追加
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
        cbar.outline.set_linewidth(0.2)
        cbar.set_ticks([toml_i["min_value"], toml_i["max_value"]])
        #cbar.locator = MaxNLocator(nbins=2)  # 目盛りを5個に設定
        #cbar.update_ticks()  # 目盛りを更新
        # 目盛り線（tick marks）を非表示にする
        cbar.ax.tick_params(axis='both', which='both', length=0)
    #print(face_colors_rgba.shape)
    #print(faces.shape)
    #print(verts.shape)
    #print(verts[faces].shape)

#    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
#        # **重心を求める**
#        centroid = np.mean(verts_transformed, axis=0)
        
        # **全頂点を重心の逆方向に移動**
        verts_transformed = verts_transformed - centroid + np.array(data_toml["view"]["translate"])


    poly3d = Poly3DCollection(verts_transformed[faces2], facecolors=face_colors_rgba, linewidths=0.01, edgecolors=face_colors_rgba, alpha=toml_i["alpha"], shade=toml_i["shade"], zorder=10)
    ax.add_collection3d(poly3d)

    #R, L_QR = project_and_distance(verts_transformed[faces2], normals[faces2], z_plane=3.0)
    #plot_heatmap(R, L_QR, face_colors_rgba, alpha=1, cmap=cmap, resolution=200, method='cubic')

    face_color_rgb = []
    for tmp in np.squeeze(face_colors_rgba):
        face_color_rgb.append(f"rgb({255*tmp[0]}, {255*tmp[1]}, {255*tmp[2]})")


    figp.add_trace(go.Mesh3d(
        x=verts_transformed[:,0],  # x座標.astype(np.float16)
        y=verts_transformed[:,1],  # y座標
        z=verts_transformed[:,2],  # z座標
        i=faces2[:,0],     # 面の1番目の頂点インデックス
        j=faces2[:,1],     # 面の2番目の頂点インデックス
        k=faces2[:,2],     # 面の3番目の頂点インデックス
        #color = face_color_rgb[index],
        facecolor = face_color_rgb[:],
        #intensity=color_values,  # 各メッシュごとの色
        #colorscale="Viridis",
        opacity=toml_i["alpha"],        # 透明度
        flatshading=False,  # フラットシェーディングを無効にする
        lighting=dict(
            ambient=0.7,   # 環境光（全体の明るさ）
            diffuse=0.9,   # 拡散光（物体の明るさ）
            specular=1.0,  # 鏡面反射（光の反射）
            roughness=0.2,  # 表面の粗さ
            fresnel=0.2     # フレネル効果
        ),
    ))

    plotly_colorscale = matplotlib_cmap_to_plotly_colorscale(cmap)
    dummy_intensity = intensity = np.linspace(toml_i["min_value"], toml_i["max_value"], 3)
    figp.add_trace(go.Mesh3d(
        x=[0, 0, 0], y=[0, 0, 0], z=[0, 0, 0],
        i=[0], j=[0], k=[0],
        intensity=dummy_intensity,
        colorscale=plotly_colorscale,
        showscale=True,
        opacity=0.0,
        colorbar=dict(
            title="",
            len=0.75,
            x=1.1,
            tickvals=[toml_i["min_value"], toml_i["max_value"]],
            ticktext=[str(toml_i["min_value"]), str(toml_i["max_value"])]
        )
    ))

    #rendering_by_blender(data_toml,verts_transformed, faces2)

    return

def compute_centroids_and_normals(verts):
    """
    Compute centroid and face normal for each triangle.

    Parameters
    ----------
    verts : (F, 3, 3) array
        Vertices of each triangle.

    Returns
    -------
    centroids : (F, 3) array
        Centroid of each triangle.
    face_normals : (F, 3) array
        Unit normal vector of each triangle.
    """
    # Centroids: mean of the 3 vertices
    centroids = verts.mean(axis=1)

    # Edge vectors
    v1 = verts[:, 1] - verts[:, 0]  # vector from vertex0->vertex1
    v2 = verts[:, 2] - verts[:, 0]  # vector from vertex0->vertex2

    # Cross product to get face normal
    face_normals = np.cross(v1, v2)

    # Normalize face normals
    norm = np.linalg.norm(face_normals, axis=1, keepdims=True)
    face_normals /= norm

    return centroids, face_normals

def project_and_distance(vertices, normals, z_plane):
    """
    Project points along normals to a plane at z=z_plane and compute distances.

    Parameters
    ----------
    vertices : (N,3) array
    normals : (N,3) array
    z_plane : float

    Returns
    -------
    R : (N,3) array
        Intersection points on plane.
    L_QR : (N,) array
        Distance from Q to R.
    """

    centroids, face_normals = compute_centroids_and_normals(vertices)

    Q = centroids
    n = face_normals
    Qz = Q[:, 2]
    nz = n[:, 2]

    mask = np.abs(nz) > 1e-12
    t = np.full(Q.shape[0], np.nan)
    t[mask] = (z_plane - Qz[mask]) / nz[mask]

    R = Q.copy()
    #R[mask] = Q[mask] + (t[mask, None] * n[mask])

    # Distance between Q and R
    L_QR = np.linalg.norm(R - Q, axis=1)

    return R, L_QR

# データ読み込み例（テスト用に乱数）
# vertices = np.random.randn(100,3)
# normals = np.random.randn(100,3)
# R, L_QR = project_and_distance(vertices, normals, z_plane=0.0)

def plot_heatmap(R, L_QR, face_colors_rgba, alpha, cmap, resolution=200, method='cubic'):
    """
    Create a heatmap of distances on the projected plane.

    Parameters
    ----------
    R : (N,3) array
        Intersection points on plane.
    L_QR : (N,) array
        Distance Q->R.
    alpha : float
        Decay rate for the exponential weight.
    resolution : int
        Number of grid points along each axis.
    method : str
        Interpolation method: 'linear', 'cubic', or 'nearest'.
    """
    R_x, R_y, R_z = R[:, 0], R[:, 1], R[:, 2]

    # Define grid
    #xi = np.linspace(R_x.min(), R_x.max(), resolution)
    #yi = np.linspace(R_y.min(), R_y.max(), resolution)
    xi = np.linspace(-5.0, 5.0, resolution)
    yi = np.linspace(-5.0, 5.0, resolution)
    XI, YI = np.meshgrid(xi, yi)

    # Compute weights
    min_val = np.min(L_QR)
    w = np.exp(-alpha * (L_QR-min_val))
    #face_colors_weighted = face_colors_rgba * w[:, None]

    # Interpolate L_QR onto grid
    ZI = griddata((R_x, R_y), R_z, (XI, YI), method=method)

    # Plot heatmap
    plt.figure(figsize=(6,6))
    #cmap = 'inferno'
    im = plt.imshow(
        #ZI, origin='lower', extent=(R_x.min(), R_x.max(), R_y.min(), R_y.max()),
        ZI, origin='lower', extent=(-5.0, 5.0, -5.0, 5.0),
        cmap='inferno', vmin=0.0, vmax=1.0, 
        aspect='equal'
    )
    plt.colorbar(im, label='Distance Q->R')
    plt.xlabel('X on plane')
    plt.ylabel('Y on plane')
    plt.title('Heatmap of Distance Q->R projected on plane')
    plt.tight_layout()
    plt.savefig("2dplot.pdf")  # save as PDF here
    plt.close()

# 使用例:
# R, L_QR = project_and_distance(vertices, normals, z_plane=0.0)
# plot_heatmap(R, L_QR)


def plot_lattice_box(figp, ax, data_toml, data_cube):

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

    box_vertices = box_vertices + np.array(origin) - centroid + np.array(data_toml["view"]["translate"]) 
    
    # ボックスの12本のエッジを定義
    edges = [
        [0, 1], [0, 2], [0, 3],  # 原点からのエッジ
        [1, 4], [1, 5],          # v1 からのエッジ
        [2, 4], [2, 6],          # v2 からのエッジ
        [3, 5], [3, 6],          # v3 からのエッジ
        [4, 7], [5, 7], [6, 7]   # 最後の3つのエッジ
    ]

    lattice, box_vertices, edges, edge_colors = create_supercell(lattice, box_vertices, faces=np.array(edges), supercell_size=data_toml["view"]["supercell"])

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

    lines = np.vstack([lines_front, lines_behind])
    x, y, z = [], [], []
    for line in lines:
        for point in line:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
        x.append(None)  # None を入れて線を区切る
        y.append(None)
        z.append(None)

    figp.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(color='black', width=1)
    ))

def matplotlib_cmap_to_plotly_colorscale(cmap, n=256):
    """
    Convert a matplotlib colormap to a Plotly colorscale.
    
    Parameters:
        cmap: A matplotlib colormap instance
        n: Number of discrete colors to sample
    
    Returns:
        List of [position, color] for Plotly colorscale
    """
    colorscale = []
    for i in range(n):
        frac = i / (n - 1)
        rgba = cmap(frac)
        rgb = tuple(int(255 * c) for c in rgba[:3])  # drop alpha
        colorscale.append([frac, f'rgb{rgb}'])
    return colorscale
