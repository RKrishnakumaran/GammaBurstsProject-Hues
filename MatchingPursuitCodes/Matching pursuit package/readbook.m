function [batoms, header]=readbook(filename, offset)
%function [batoms, header]=readbook(filename, offset)

if offset<0
   error('offset <= 0');
end;

SCALE  =1;
FREQ   =2;
POS    =3;
MODULUS=4;
AMPLI  =5;

H_SAMPLING_FREQ=1;
H_SIGNAL_SIZE=2;
H_POINTS_PER_MICROVOLT=3;
H_VERSION=4;

if ~exist(filename, 'file')
   error([filename ' does not exist']);
end

[batoms ,size, sampling, file_offset, header, book_version]=readrawb(filename, offset);

batoms(:,SCALE)   = batoms(:,SCALE)./header(H_SAMPLING_FREQ);
batoms(:,FREQ)    = batoms(:,FREQ).*header(H_SAMPLING_FREQ)/header(H_SIGNAL_SIZE);
batoms(:,POS)     = batoms(:,POS)./header(H_SAMPLING_FREQ);
batoms(:,MODULUS) = batoms(:,MODULUS)./header(H_POINTS_PER_MICROVOLT);
batoms(:,AMPLI)   = 2.*batoms(:,AMPLI)./header(H_POINTS_PER_MICROVOLT);
