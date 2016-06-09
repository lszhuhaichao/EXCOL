clear all
close all

addpath('GCO/');
addpath('shape/');

%%%%%%%%%%%%%%shape context parameters%%%%%%%%%%%%%%
global n_dist;
global n_theta;
global bTangent;
global bSmoothCont;
global n_contsamp;

n_dist		= 5;
n_theta		= 36;
bTangent	= 1;
bSmoothCont	= 100;
n_contsamp	= 100;

%%%%%%%%%%%%%% EXCOL parameters %%%%%%%%%%%%%%%%%%%%
global lambda1;
global lambda2;
global lambda3;

lambda1  = 0.00001;
lambda2 = 100;
lambda3 = 1;

%%%%%%%%%%%%%% case parameters %%%%%%%%%%%%%%%%%%%%
folder = 'mouse/';
frames = 6:18;

%%%%%%%%%%%%%% EXCOL Begin %%%%%%%%%%%%%%%%%%%%
load(['examples/' folder 'patches_o/' int2str(frames(1)) '.mat']);
patches_o = patches;
load(['examples/' folder 'patches/' int2str(frames(1)) '.mat']);
load(['examples/' folder 'descriptors/' int2str(frames(1)) '.mat']);
im1 = imread(['examples/' folder 'original_frame/' int2str(frames(1)) '.png']);
patches1 = patches;
localDescriptor1 = localDescriptor;

figure();
imshow(im1);
figure();
imshow(patches_o,[]);
if(~exist(['examples/' folder 'excol/init' int2str(frames(1)) '.mat'],'file'))
    mkdir(['examples/' folder 'excol/']);
    NumLabels = input('Input Number of Layers?');
    regionsLabels1 = NumLabels*ones(size(unique(patches)));
    while(1)
        for i=1:NumLabels-1
            disp(['Please set ' int2str(i) 'th Layer.']);
            [x,y] = ginput();
            for cc = 1:length(x)
                regionsLabels1(patches(round(y(cc)),round(x(cc)))) = i;
            end
        end
        layerMap = zeros(size(patches));
        for i=1:length(regionsLabels1)
            layerMap(patches_o==i)=regionsLabels1(i);
        end
        figure();
        imshow(label2rgb(layerMap));
        isCorrect = input('The Layer is Correct?(true/false)');
        if(isCorrect)
            break;
        end
    end
    save(['examples/' folder 'excol/init' int2str(frames(1)) '.mat'],'regionsLabels1');
else
    load(['examples/' folder 'excol/init' int2str(frames(1)) '.mat']);
    NumLabels = length(unique(regionsLabels1));
    layerMap = zeros(size(patches));
    for i=1:length(regionsLabels1)
        layerMap(patches_o==i)=regionsLabels1(i);
    end
    figure();
    imshow(label_to_color(layerMap));
end

mkdir(['examples/' folder 'excol/result/']);
imwrite(label_to_color(layerMap),['examples/' folder 'excol/result/' int2str(frames(1)) '.png']);
figure();
for f=2:length(frames)
    disp(['Propagate the layer labels from frame ' int2str(frames(1)) ' to ' int2str(frames(f))]);
    load(['examples/' folder 'patches_o/' int2str(frames(f)) '.mat']);
    patches_o = patches;
    load(['examples/' folder 'patches/' int2str(frames(f)) '.mat']);
    load(['examples/' folder 'descriptors/' int2str(frames(f)) '.mat']);
    load(['examples/' folder 'adjacency/' int2str(frames(f)) '.mat']);
    
    NumSites = length(unique(patches));
    h = GCO_Create(NumSites, NumLabels);
    
%%%%%%%%%  Set The DataCost %%%%%%%%%%%%%%%%%%%%%%%    
    disp('Set The Data Cost ...');
    DataCost = zeros(NumLabels,NumSites);
    for i=1:NumSites
        for j=1:NumLabels
            inds = find(regionsLabels1==j);
            subDeiscriptors = localDescriptor1(inds);
            
            ms = zeros(size(inds));
            for jj=1:length(inds)
                m = mSimilarity(localDescriptor1{inds(jj)}, localDescriptor{i});
                ms(jj) = m;
            end
            [vv,ii] = min(ms);
            DataCost(j,i) = vv;
        end
    end
    DataCost = int32(round(lambda2*bsxfun(@rdivide, DataCost, sum(DataCost,1))));
    GCO_SetDataCost(h,DataCost);
    
    disp('Set The Smooth Cost ...');
    SmoothCost = 1-eye(NumLabels);
    GCO_SetSmoothCost(h,SmoothCost);
    
    disp('Set The Neighbors ...');
    Neighbors = zeros(NumSites);
    for j=1:size(adjPairs,1)
        Neighbors(adjPairs(j,1),adjPairs(j,2)) = 1/(1+norm(localDescriptor{adjPairs(j,1)}.colorHist-localDescriptor{adjPairs(j,2)}.colorHist));
    end
    Neighbors = round(lambda3*Neighbors);
    GCO_SetNeighbors(h, Neighbors);
    
    disp('Graph Cut Optimization Swap ...');
    GCO_Swap(h,100);
    regionsLabels = GCO_GetLabeling(h);
    
    
    layerMap = zeros(size(patches));
    for i=1:length(unique(patches))
        layerMap(patches_o==i)=regionsLabels(i);
    end
    
    imshow(label_to_color(layerMap));
    imwrite(label_to_color(layerMap),['examples/' folder 'excol/result/' int2str(frames(f)) '.png']);
    save(['examples/' folder 'excol/result/' int2str(frames(f)) '.mat'],'regionsLabels');
end