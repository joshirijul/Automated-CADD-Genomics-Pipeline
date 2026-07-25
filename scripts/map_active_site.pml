# PyMol Macro: BACE1 Active Site Topography & Pocket Mapping
# Usage in terminal: pymol -c map_active_site.pml (headless) OR pymol map_active_site.pml (GUI)

reinitialize
load 1fkn.pdb, bace1

# Remove water molecules to clean up the workspace
remove resn HOH

# Set background to clean white for report-ready presentation
bg_color white
hide all

# Show secondary structure cartoon for overall enzyme
show cartoon, bace1
color slate, bace1
set cartoon_transparency, 0.2, bace1

# Isolate the co-crystallized ligand (heteroatoms excluding water)
select bound_ligand, hetatm and bace1
show sticks, bound_ligand
util.cbay bound_ligand

# MAP THE ACTIVE SITE POCKET: Select residues within 5 Angstroms of the ligand
select active_site_pocket, byres (bound_ligand around 5.0)
show sticks, active_site_pocket
color lightpink, active_site_pocket

# Render transparent surface topography specifically over the binding pocket
show surface, active_site_pocket
set surface_color, white, active_site_pocket
set transparency, 0.4, active_site_pocket

# Orient view to focus on the active site and cast high-def shadows
zoom active_site_pocket, 8.0
set ray_shadows, 1
set orthoscopic, on

# Render high-resolution PNG for git repo and lab report
ray 1600, 1200
png bace1_active_site_topography.png
quit