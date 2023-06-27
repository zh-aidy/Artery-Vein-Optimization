# File Structure
```
Artery-Vein-Optimization
├─ data
│  ├─ AV-DRIVE
│  │  ├─ optim_result
│  │  │  └─ image_vessel_correct.png
│  │  ├─ orig
│  │  │  ├─ 1_image.png
│  │  │  └─ 1_mask.png
│  │  └─ pred_result_cnn
│  │     ├─ 1_pred_label.png
│  │     └─ 1_pred_vessel.png
│  ├─ CY-BIT
│  │  ├─ optim_result
│  │  │  └─ image_vessel_correct.png
│  │  ├─ orig
│  │  │  ├─ 5126_image.png
│  │  │  ├─ 5126_mask.png
│  │  │  └─ CY-BIT.zip
│  │  └─ pred_result_cnn
│  │     ├─ 5126_pred_label.png
│  │     └─ 5126_pred_vessel.png
│  └─ LES
│     ├─ optim_result
│     │  └─ image_vessel_correct.png
│     ├─ orig
│     │  ├─ 3_image.png
│     │  └─ 3_mask.png
│     └─ pred_result_cnn
│        ├─ 3_pred_label.png
│        └─ 3_pred_vessel.png
├─ main.py
├─ README.md
└─ src
   ├─ color_diff.py
   ├─ function_label_correct.py
   ├─ function_pruning.py
   ├─ function_topolpgy_adjustment.py
   ├─ function_trace_source.py
   ├─ MulSke_graph_optimize.py
   ├─ mul_ske.py
   ├─ Nearest_neighbor_label_passing.py
   └─ shorted_path_dijkstra.py
```
## Introduction

Article title:

**Optimization of retinal artery/vein classification based on vascular topology**

The key idea of our proposed post-processing optimization method is using the prior knowledge of tree-structure artery and vein, for example, only one label within one consecutive vascular segment, the label of vessel branch and the one vessel trunk should be same, and the blood flow direction is from the optic disc to the vessel end.

The CY-BIT dataset annotated by ophthalmologists from Beijing Chaoyang Hospital, contains 56 retinal images from Chinese patients with different disease. These retinal images contain a variety of clinical conditions that affect the clarity of retinal vessels, such as uneven illumination, clear choroid, hemorrhages, exudations and so on. The training set comprises 38 retinal images with $30^{\circ}$ FOV and $1444\times 1620$ pixels (29 images), and $45^{\circ}$ FOV and $1958\times 2196$ pixels (9 images). The test set comprises 18 retinal images with $30^{\circ}$ FOV and $1444\times 1620$ pixels (16 images), and $45^{\circ}$ FOV and $1958\times 2196$ pixels (2 images).

We provide one sample image for three datasets respectively. 

The prediction results in `pred_result_cnn` folder are obtained by training the VCNet[1], if you waant to test other images, you can download the `.pkl` file though following https URLs:

* AV-DRIVE
[https://pan.baidu.com/s/1HuYERSPr0n0zxA1b3T2lgQ](https://pan.baidu.com/s/1HuYERSPr0n0zxA1b3T2lgQ 
) password: h30i

* LES
[https://pan.baidu.com/s/1fn4n9Zfe08Z2lhIIZU0JSg](https://pan.baidu.com/s/1fn4n9Zfe08Z2lhIIZU0JSg 
) password: tljf

* CY-BIT
[https://pan.baidu.com/s/1pfE6ykJafO18o-F9RfLHOA](https://pan.baidu.com/s/1pfE6ykJafO18o-F9RfLHOA 
) password: hlk4

# Run
You can run the demo according to following code:
* AV-DRIVE

```
python main.py --fundus_img_path data\AV-DRIVE\orig\1_image.png --pred_av_path data\AV-DRIVE\pred_result_cnn\1_pred_label.png --pred_v_path data\AV-DRIVE\pred_result_cnn\1_pred_vessel.png --opti_save_path data\AV-DRIVE\optim_result --H (290, 485) --Ks 1.
```

* LES

```
python main.py --fundus_img_path data\LES\orig\3_image.png --pred_av_path data\LES\pred_result_cnn\3_pred_label.png --pred_v_path data\LES\pred_result_cnn\3_pred_vessel.png --opti_save_path data\LES\optim_result --H (720, 803) --Ks 2.5
```

* CY-BIT

```
python main.py --fundus_img_path data\CY-BIT\orig\5126_image.png --pred_av_path data\CY-BIT\pred_result_cnn\5126_pred_label.png --pred_v_path data\CY-BIT\pred_result_cnn\5126_pred_vessel.png --opti_save_path data\CY-BIT\optim_result --H (745, 751) --Ks 2.5
```

##  Annotations
The parameter $H$ is the center point coordinates of the optic disc, so $H$ is different in each fundus image. We get this parameter using [another project](https://doi.org/10.1016/j.compbiomed.2023.106796) of our team. 

# Reference
```bibtex
@article{sh20,
  title={Automatic artery/vein classification using a vessel-constraint network for multicenter fundus images},
  author={Hu, Jingfei and Wang, Hua and Cao, Zhaohui and Wu, Guang and Jonas, Jost B and Wang, Ya Xing and Zhang, Jicong},
  journal={Frontiers in cell and developmental biology},
  pages={1194},
  year={2021},
  publisher={Frontiers}
}
```
