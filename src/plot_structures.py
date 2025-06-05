import os
from pymol import cmd as pymol_cmd
from pymol import util


def plot_structure(predicted_structure_superposed_filepath, output_filepath, orient_structure=True):
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
    pymol_cmd.ray(width=2400, height=2400)
    pymol_cmd.png(output_filepath, dpi=300)


def plot_mmalign_superposed_structures(predicted_structure_superposed_filepath, native_structure_filepath,
                                       output_filepath, orient_structure=True):
    pymol_cmd.reinitialize()
    pymol_cmd.load(predicted_structure_superposed_filepath)
    pymol_cmd.load(native_structure_filepath)
    obj_list = pymol_cmd.get_names('objects')
    util.cbc(obj_list[0]) # color predicted structure by chain
    pymol_cmd.color('grey80', obj_list[1]) # native structure in gray
    pymol_cmd.set('cartoon_transparency', 0.2, obj_list[1])
    pymol_cmd.show("cartoon")
    pymol_cmd.hide("lines", "all")
    pymol_cmd.hide("spheres", "all")
    pymol_cmd.hide('sticks', "all")
    pymol_cmd.hide('nonbonded', "all")
    pymol_cmd.set("depth_cue", 0)
    pymol_cmd.set("spec_reflect", 0)
    pymol_cmd.space('cmyk')
    pymol_cmd.set('opaque_background', 'on')
    pymol_cmd.set('bg_rgb', [255, 255, 255])
    if orient_structure:
        pymol_cmd.orient()
    pymol_cmd.zoom(complete=1)
    pymol_cmd.ray(width=2400, height=2400)
    session_output = os.path.splitext(output_filepath)[0] + '.pse'
    pymol_cmd.save(session_output)
    pymol_cmd.png(output_filepath, dpi=300)


def plot_from_session(pymol_session_filepath, output_filepath, orient_structures=False):
    pymol_cmd.reinitialize()
    pymol_cmd.load(pymol_session_filepath)
    pymol_cmd.set('bg_rgb', [255, 255, 255])
    if orient_structures:
        pymol_cmd.orient()
    pymol_cmd.zoom(complete=1)
    pymol_cmd.ray(width=2400, height=2400)
    pymol_cmd.png(output_filepath, dpi=300)


def plot_from_session_highlight_interface(pymol_session_filepath, output_filepath, pathogen_chain, pathogen_if_residues,
                                          human_interactor_chain, human_interactor_if_residues,
                                          orient_structures=False):
    pymol_cmd.reinitialize()
    pymol_cmd.load(pymol_session_filepath)
    obj_list = pymol_cmd.get_names('objects')
    pymol_cmd.set('bg_rgb', [255, 255, 255])
    pymol_cmd.set('opaque_background', 'on')
    pymol_cmd.set('orthoscopic', 'on')
    pymol_cmd.set('cartoon_transparency', 0.65)
    pathogen_if_residues_str = '+'.join([str(x) for x in pathogen_if_residues])
    human_interactor_if_residues_str = '+'.join([str(x) for x in human_interactor_if_residues])
    pymol_cmd.set('cartoon_transparency', 0.0, f'{obj_list[0]} & chain {pathogen_chain} & resi {pathogen_if_residues_str}')
    pymol_cmd.set('cartoon_transparency', 0.0, f'{obj_list[1]} & chain {human_interactor_chain} & resi {human_interactor_if_residues_str}')
    if orient_structures:
        pymol_cmd.orient()
    pymol_cmd.ray(width=2400, height=2400)
    pymol_cmd.png(output_filepath, dpi=300)

