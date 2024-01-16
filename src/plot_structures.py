from pymol import cmd as pymol_cmd
from pymol import util
from pathlib import Path


def plot_structure(pdb_path: Path, out_path: Path) -> None:
    """Display 3D structure"""
    pymol_cmd.load(str(pdb_path))
    pymol_cmd.spectrum("b", "red_white_blue", "*", 0, 100)
    pymol_cmd.show("cartoon")
    pymol_cmd.hide("lines", "all")
    pymol_cmd.set("depth_cue", 0)
    pymol_cmd.set("spec_reflect", 0)
    pymol_cmd.set("bg_rgb", [1, 1, 1])
    pymol_cmd.set("orthoscopic", "on")
    pymol_cmd.set("cartoon_fancy_helices", 1)
    pymol_cmd.set("cartoon_smooth_loops", 1)
    pymol_cmd.set("cartoon_highlight_color", 1)
    pymol_cmd.set("cartoon_loop_radius", 0.5)
    pymol_cmd.set("ray_shadows", 0)
    pymol_cmd.set("ray_texture", 1)
    pymol_cmd.ray()
    pymol_cmd.png(str(out_path))


def plot_mmalign_superposed_structures(predicted_structure_superposed_filepath, native_structure_filepath,
                                       output_filepath, orient_structure=True):
    pymol_cmd.reinitialize()
    pymol_cmd.load(predicted_structure_superposed_filepath)
    pymol_cmd.load(native_structure_filepath)
    obj_list = pymol_cmd.get_names('objects')
    util.cbc(obj_list[0]) # color predicted structure by chain
    pymol_cmd.color('grey80', obj_list[1]) # native structure in gray
    pymol_cmd.show("cartoon")
    pymol_cmd.hide("lines", "all")
    pymol_cmd.hide("spheres", "all")
    pymol_cmd.hide('sticks', "all")
    pymol_cmd.hide('nonbonded', "all") # https://pymolwiki.org/index.php/Hide
    pymol_cmd.set("depth_cue", 0)
    pymol_cmd.set("spec_reflect", 0)
    pymol_cmd.space('cmyk')
    pymol_cmd.set('opaque_background', 'on')
    pymol_cmd.set('bg_rgb', [255, 255, 255])
    if orient_structure:
        pymol_cmd.orient()
    pymol_cmd.zoom(complete=1)
    #pymol_cmd.move(string axis, float distance)
    #pymol_cmd.turn(string axis, float angle)
    #pymol_cmd.rotate(list-or-string axis, angle=0, string selection = "all", int state = 0, int camera = 1, string object = None)
    pymol_cmd.ray(width=2400, height=2400)
    pymol_cmd.png(output_filepath, dpi=300)


def plot_struct(predicted_structure_superposed_filepath, output_filepath, orient_structure=True):
    pymol_cmd.reinitialize()
    pymol_cmd.load(predicted_structure_superposed_filepath)
    obj_list = pymol_cmd.get_names('objects')
    util.cbc(obj_list[0])  # color predicted structure by chain
    pymol_cmd.show("cartoon")
    pymol_cmd.hide("lines", "all")
    pymol_cmd.hide("spheres", "all")
    pymol_cmd.hide('sticks', "all")
    pymol_cmd.hide('nonbonded', "all")  # https://pymolwiki.org/index.php/Hide
    pymol_cmd.set("depth_cue", 0)
    pymol_cmd.set("spec_reflect", 0)
    pymol_cmd.space('cmyk')
    pymol_cmd.set('opaque_background', 'on')
    pymol_cmd.set('bg_rgb', [255, 255, 255])
    if orient_structure:
        pymol_cmd.orient()
    pymol_cmd.zoom(complete=1)
    # pymol_cmd.move(string axis, float distance)
    # pymol_cmd.turn(string axis, float angle)
    # pymol_cmd.rotate(list-or-string axis, angle=0, string selection = "all", int state = 0, int camera = 1, string object = None)
    pymol_cmd.ray(width=2400, height=2400)
    pymol_cmd.png(output_filepath, dpi=300)
