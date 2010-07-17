
gfx r n Region_1
gfx r e Region_1
gfx def faces egroup Region1

gfx change_identifier node +100000
gfx change_identifier ele +100000

gfx r n Region_2
gfx r e Region_2
gfx def faces egroup Region2

gfx change_identifier node +100000
gfx change_identifier ele +100000

gfx r n Interface
gfx r e Interface
gfx def faces egroup Interface_set
