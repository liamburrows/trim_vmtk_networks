This code runs in Python within the 3D Slicer environment, accessed via Jupyter Notebook using the SlicerJupyter extension.


### Setup
1. Install 3D Slicer and the above extensions via the Extension Manager.
2. Launch a Jupyter server from within Slicer (via SlicerJupyter).
3. Connect to the Slicer kernel in Jupyter Notebook — **do not use a standard Python kernel**.
4. Clone this repository and open `extract_trim_network.ipynb`.

## Usage

1. Load your segmented vascular label map (e.g. `example_patree.nii.gz`) into the Slicer scene.
2. Open and run `extract_trim_network.ipynb`, specifying:
   - Segmented vessel file_ename (`example_patree.nii.gz`)
   - The maximum number of bifurcations (`N`)
3. The mesh, network, root and trimmed network is saved locally.


## Example

<table>
  <tr>
    <td align="center"><img src="ims/mesh.png" width="400"/><br/><b>Figure 1:</b> Mesh</td>
    <td align="center"><img src="ims/net.png" width="400"/><br/><b>Figure 2:</b> Network from vmtk</td>
  </tr>
  <tr>
    <td align="center"><img src="ims/root.png" width="400"/><br/><b>Figure 3:</b> Root of artery</td>
    <td align="center"><img src="ims/trim.png" width="400"/><br/><b>Figure 4:</b> Trimmed network</td>
  </tr>
</table>
