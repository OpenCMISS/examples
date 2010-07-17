
gfx r n TwoRegion1.part0.exnode
#gfx r e TwoRegion1.part0.exelem generate

gfx change_identifier node_offset +1000
gfx change_identifier element_offset +1000

gfx r n TwoRegion2.part0.exnode
#gfx r e TwoRegion2.part0.exelem generate

gfx change_identifier node_offset +1000
gfx change_identifier element_offset +1000

gfx r n TwoRegion3.part0.exnode
#gfx r e TwoRegion3.part0.exelem generate
