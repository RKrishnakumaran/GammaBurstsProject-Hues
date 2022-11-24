% first create signal with some Gabor functions, using 'gabor.m' 
% called with the following parameters:  
% gabor(_size,_sampling, width, frequency, position, amplitude, phase)
l=4;
f=128;
s=zeros(1,l*f);
s=s+gabor(l,f, 0.05, 1,  .3, 5, -.3);
s=s+gabor(l,f, 0.25,30,   1, 1,   0);
s=s+gabor(l,f, 0,    0, 2.5, 2,   0);
s=s+gabor(l,f, 0.25,10, 3.5, 1,   0);
s=s+gabor(l,f, 4,   40, 2,  .1,   0);


figure(1);subplot(211);
t=(1:l*f)./f;
plot(t,s); title('original signal');

%save the signal to a file
save tmp.txt s -ascii

%and decomposed using mp4 with commands:
%MP4>set -e 95 -i 10
%MP4>reinit -O 512 -R 1000000
%MP4>loadsig -O tmp.txt -F 128
%MP4>mp
%MP4>save

% now we load the results of the decomposition, atoms in 'a' and header in 'h':

[a, h]=readbook('book.b', 0);

% some constants related to the internal structure of the 'book':
SCALE  =1;
FREQ   =2;
POS    =3;
MODULUS=4;
AMPLI  =5;
PHASE  =6;
H_SAMPLING_FREQ=1;
H_SIGNAL_SIZE=2;
H_POINTS_PER_MICROVOLT=3;
H_VERSION=4;

% .. so we can reconstruct the signal from the 'atoms' stored in the book

rec=zeros(1,h(H_SIGNAL_SIZE));

for i=1:size(a,1)
rec=rec+...
gabor(h(H_SIGNAL_SIZE)/h(H_SAMPLING_FREQ),h(H_SAMPLING_FREQ),a(i,SCALE),a(i,FREQ),a(i,POS),a(i,AMPLI),a(i,PHASE));end

t1=(1:h(H_SIGNAL_SIZE))./h(H_SAMPLING_FREQ);
figure(1); subplot(212); plot(t1,rec); 
title(sprintf('reconstruction from %d MP iterations', size(a,1)));

% ... and plot the t-f representation

[map,xx,yy]=mp2tf(a, h, 1, 1);
figure(2);
imagesc(xx,yy,map); set(gca,'ydir', 'normal')
