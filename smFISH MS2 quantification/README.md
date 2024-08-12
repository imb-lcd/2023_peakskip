Run `procedure6_pxprob_curate.m` to quantify smFISH and MS2 signal intensity.

* `procedure6_pxprob_curate.m` calls `get_value_seg.m` to get the value of segmented nuclei from the data in the `Track data` folder.
* `procedure6_pxprob_curate.m` calls `pixel_filter5_pxprob_curate.m` to get the smFISH and MS2 signal (`ms2_array` and `smFISH_array`).

The result is saved in the `result` folder.
