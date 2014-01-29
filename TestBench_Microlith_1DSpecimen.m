%% Test computation of 1D specimen with 2D system.

% Compare images computed by 1D sum over source and 2D sum over source.

% Specimen: double-slit.
x=[-3:0.01:3];
slitwidth=0.3;
specimen1d= single(abs(x-0.5*slitwidth)<1E-10 | abs(x+0.5*slitwidth)<1E-10);
specimen2d= single(repmat(specimen1d,length(x),1));

%% Compute 2D and 1D image of fluorescence microscope.
params.NAo=0.95;
params.NAc=0.9;
params.wavelength=0.55;
params.nEmbb=1;
params.nImm=1;



fluor=microlith(x,0);
fluor.computesys('Fluorescence',params);
fluor.computeimage(specimen2d,'CPU');



fluor1d=microlith(x,0);
fluor1d.computesys('Fluorescence',params);
fluor1d.computeimage(specimen1d,'CPU');


%% Compute DIC image of phase object.

params.NAo=0.95;
params.NAc=0.9;
params.wavelength=0.55;
params.nEmbb=1;
params.nImm=1;
params.shearangle=0;
params.shear=0.2;
params.bias=pi/4;

tic2D=tic;
DIC2d=microlith(x,0);
DIC2d.computesys('DIC',params);
DIC2d.Ic(DIC2d.mm>0)=0;
DIC2d.computeimage(exp(1i*specimen2d),'CPU');
time2D=toc(tic2D);

tic1D=tic;

DIC1d=microlith(x,0);
DIC1d.computesys('DIC',params);
DIC1d.Ic(DIC1d.mm>0)=0;
DIC1d.computeimage(exp(1i*specimen1d),'CPU');

time1D=toc(tic1D);


%% Plot the results.



figure(1); clf; set(1,'Position',[100 100 800 1200],'defaultaxesfontsize',12,'color','white');

subplot(321); plot(x,specimen1d); title('Transmission'); xlim([-1 1]);
subplot(322); imagesc(x,x,specimen2d); axis equal; title('Transmission in 2D'); xlim([-1 1]);

subplot(324); imagesc(x,x,fluor.img); axis equal; title('Fluorescence image');  xlim([-1 1]);
subplot(323); plot(x,gray2norm(fluor.img(1,:)),':',x,gray2norm(fluor1d.img),'--','LineWidth',2); 
legend('Profile from 2D image','1D image', 'Location', 'NorthOutside');


subplot(326); imagesc(x,x,DIC2d.img); axis equal; title('DIC image');  xlim([-1 1]);
subplot(325);
plot(x,gray2norm(DIC2d.img(1,:)),'-.',x,gray2norm(DIC1d.img),'--','LineWidth',2);
legend(['Profile through 2D image: ' num2str(time2D) 's'],...
    ['Directly computed 1D profile:' num2str(time1D) 's'],'Location','NorthOutside'); 