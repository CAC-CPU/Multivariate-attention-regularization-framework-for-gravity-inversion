function h=corob(x,y,z,cs,col)

% Patch 3-D block for the model

x0=x-cs(1)/2;
x1=x+cs(1)/2;
y0=y-cs(2)/2;
y1=y+cs(2)/2;
z0=z-cs(3)/2;
z1=z+cs(3)/2;
Ver=[x0 y0 z0;
    x1 y0 z0;
    x1 y1 z0;
    x0 y1 z0;
    x0 y0 z1;
    x1 y0 z1;
    x1 y1 z1;
    x0 y1 z1];

Fac=[1 2 6 5;
    2 3 7 6;
    3 4 8 7;
    4 1 5 8;
    1 2 3 4;
    5 6 7 8];
h=patch('Vertices',Ver,'Faces',Fac,'Edgecolor',[0.4 0.4 0.4],'FaceColor',col);