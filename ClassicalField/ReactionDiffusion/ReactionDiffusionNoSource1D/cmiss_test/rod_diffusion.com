fem def coor 1,1
fem def base;r;linear
fem def nod;r;rod
fem def ele;r;rod


fem def equa;r;diffusion
fem def mate;r;diffusion
fem def init;r;diffusion
fem def solv;r;diffusion

fem solv
