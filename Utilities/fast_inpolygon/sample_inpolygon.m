
function []=sample_inpolygon()

gridWidth=640;
gridHeight=640;

% create a grid for testing
[X,Y]=meshgrid(1:gridWidth, 1:gridHeight);
x=X(:);
y=Y(:);

% create a polygon
numPts=80;
L = linspace(0, 2.*pi, numPts); 
xv = cos(L)';
yv = sin(L)';
xv = [xv ; xv(1)]; 
yv = [yv ; yv(1)];

% map the polygon to our grid: It will be outside the array as well!
offset=40;
xv=(xv-min(xv)).*gridWidth/3+offset;
yv=(yv-min(yv)).*gridHeight/3+offset;

% test Darren's inpoly function
tic(); [In1 On1]=inpoly([x y],[xv yv]); toc();

% test my inpolygon_fast function
% Notice that it uses a similiar call interface to inpolygon not to inpoly
tic(); [In2 On2]=inpolygon_fast(x, y, xv, yv); toc();

% display the performance comparison
disp('Error:');
% this error should exactly be 0
error = sum(abs(double(In1)-double(In2))) + sum(abs(double(On1)-double(On2)))

% display the results: Display a filled polygon
figure, plot(x(find(In2==1)), y(find(In2==1)),'r.' );
hold on,  plot(xv,yv,'b-');

end
