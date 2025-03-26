
import os
import pprint
import copy
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator

import plot_style
from data_cube import Data_cube
from plot_data import plot_structure
from plot_data import plot_isosurface
from plot_data import plot_lattice_box


def plot_figures(data_toml):
    # 3Dプロット

    setting = plot_style.Plot_style_1d()
    setting.plt_style()

    # Set the backend to TkAgg
    matplotlib.use("TkAgg")

    fig = plt.figure(figsize=(7, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])  # X:Y:Z のスケールを等しくする

    # Get the Tk window manager object
    manager = plt.get_current_fig_manager()
    wposition_x = 0
    wposition_y = 0
    manager.window.wm_geometry(f"+{wposition_x}+{wposition_y}")  # Set window position (x=100, y=100)

    # 軸ラベル
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ############  データのプロット  ##############
    data_cube1 = Data_cube()                                                                                                                                         
    data_cube1.read(data_toml["cube1"])
    data_cube2 = Data_cube()                                                                                                                                         
    data_cube2.read(data_toml["cube2"])

    if(data_toml["show_structure"]):
        plot_structure(ax, data_toml, data_cube1)

    if(data_toml["show_isosurface"]):
        plot_isosurface(ax, data_toml, data_cube1, data_cube2)

    if (data_toml["show_lattice"]):
        plot_lattice_box(ax, data_toml, data_cube1)

    ######################################


    ############  set plot range  ##############
    # 各軸の範囲を取得
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()

    x_min = x_min/data_toml["view"]["ratio_zoom"]
    y_min = y_min/data_toml["view"]["ratio_zoom"]
    z_min = z_min/data_toml["view"]["ratio_zoom"]
    x_max = x_max/data_toml["view"]["ratio_zoom"]
    y_max = y_max/data_toml["view"]["ratio_zoom"]
    z_max = z_max/data_toml["view"]["ratio_zoom"]
    
    # ★ 最大の範囲を求めて、すべての軸に適用
    #ratio_zoom = 1.5
    max_range = (max(x_max - x_min, y_max - y_min, z_max - z_min) / 2)
    #max_range = (max(x_max - x_min, y_max - y_min, z_max - z_min) / 2)/data_toml["view"]["ratio_zoom"]
    
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
        print("show plot start")
        plt.show(block=True)
        print("show plot end")
        #plt.pause(1.0)
        #plt.draw()
        #plt.pause(1.0)

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



