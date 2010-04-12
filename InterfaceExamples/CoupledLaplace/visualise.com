
gfx r n TwoRegion_1.part0.exnode
gfx r e TwoRegion_1.part0.exelem generate

gfx change_identifier node_offset +1000
gfx change_identifier element_offset +1000

gfx r n TwoRegion_2.part0.exnode
gfx r e TwoRegion_2.part0.exelem generate



