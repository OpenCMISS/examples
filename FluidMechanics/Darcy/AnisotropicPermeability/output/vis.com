gfx read node './cmgui.exnode'
gfx read elem './cmgui.exelem'
gfx create window 1
#gfx mod win 1 layout width 1024 height 1024

gfx edit scene
gfx edit spectrum

