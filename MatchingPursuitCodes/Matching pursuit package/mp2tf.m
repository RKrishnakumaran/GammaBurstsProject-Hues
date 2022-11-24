function [wigXY,xx,yy]=mp2tf(book, header, Dt, Df, minF, maxF, minT, maxT)
%function [wigXY,xx,yy]=mp2tf(book, header, Dt, Df, minF, maxF, minT, maxT)
% minF, maxF, minT, maxT, dt, df -- in samples

H_SAMPLING_FREQ=1;
H_SIGNAL_SIZE=2;

dimBase=header(H_SIGNAL_SIZE);
Fsamp=header(H_SAMPLING_FREQ);

if nargin<5
	 minF=0;
	 maxF=dimBase/2;
	 minT=0;
	 maxT=dimBase;
end
if nargin<3
         Dt=1;
	 Df=1;
end


t=1:maxT-minT; % skala czasu w punktach
f=1:maxF-minF; % skala czestosci w punktach
Tsize=length(t);
Fsize=length(f);
Xsize=(maxT-minT)/Dt;
Ysize=(maxF-minF)/Df;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CAUTION: the values below are correct only if the MP decomposition
% was calculated with border conditions as zeros outside 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
book_min_t_points = 0;
book_max_t_points = dimBase;
book_min_f_points = 0;
book_max_f_points = dimBase/2;


[wigXY, xx, yy] = deal([],[],[]);

% sprawdzamy czy wyjsciowe piksle zawieraja calkowita liczbe pikseli dokladnej mapy

if mod(Tsize,Xsize)~=0
   disp(sprintf('wrong Xsize, possible values are: \n'));
   for i=1:Tsize
      if mod(Tsize,i)==0
         disp(i)
      end
   end
   error('... so the end');
elseif mod(Fsize,Ysize)~=0
   disp(sprintf('wrong Ysize, possible values are: \n'));
   for i=1:Fsize
      if mod(Fsize,i)==0
         disp(i)
      end
   end   
   error('... so the end');
end

% tlow = (t-1); %granice calkowania - male piksele
% tupp = t;  %granice calkowania - male piksele
% 
% flow = (f-1); %granice calkowania - male piksele
% fupp = f; %granice calkowania - male piksele
% 
% wig=zeros(Tsize,Fsize);
% wigTF=zeros(Tsize,Fsize);

DX=Tsize/Xsize;
DY=Fsize/Ysize;
tx=minT+t(1:DX:Tsize); %brzegi duzych pikseli na osi czasu
fy=minF+f(1:DY:Fsize); %brzegi duzych pikseli na osi czestosci

tlowx = ((tx-1)); %granice calkowania - duze piksele
tuppx = ( tx-1 + DX); %granice calkowania - duze piksele

flowy = ((fy-1));  %granice calkowania - duze piksele
fuppy = ( fy-1 + DY);  %granice calkowania - duze piksele

xx = tx+DX/2; % srodki duzych pikseli na osi czasu (T==>X)
yy = fy+DY/2; % srodki duzych pikseli na osi czestosci (F==>Y)

df=dimBase/Fsamp;% konwersja Hz na punkty

PI4=4.0*pi;
SQRT_PI = sqrt(pi);
BI_SQRT_PI = 2*SQRT_PI;
num_atom=size(book, 1);
DX_DY = DX*DY;

% dalej = 0;

wigXY=zeros(Xsize,Ysize);
wig_gabXY=sparse([]);

hm = floor(num_atom/10);

for k=1:num_atom
   %    disp(sprintf('atom %d',k));
   modulus=book(k,4);
   amplitude=book(k,5)/2;
   scale=book(k,1)*Fsamp; % skala atomu  w punktach
   freq=(df*book(k,2));   % czestosc atomu w punktach
   trans=(book(k,3)*Fsamp); %pozycja atomu w punktach
   u=book(k,3); %u=trans/Fsamp;
   
   if scale~=0
      if scale==dimBase % sinus
         %          freqTF=1+floor(freq-minF);  % przesuwamy czestosc - floor -zeby pasowala do siatki wig, +1 bo w matlabie maceirze maja indeksy od 1
         freqXY=fy(1+floor((freq-minF)/DY));  % przesuwamy czestosc - floor -zeby pasowala do siatki wig, +1 bo w matlabie maceirze maja indeksy od 1
         if (freqXY>=1) & (freqXY<=Ysize)
            %             wigTF(:,freqTF)=wigTF(:,freqTF)+(modulus^2)/Tsize; % normalizacja HAK
            wigXY(:,freqXY)=wigXY(:,freqXY)+(modulus^2)/Tsize; % normalizacja HAK
         end;
      else %gabor  
         %f_scale=PI4*(dy*scale/(2.0*pi))^2;
         %             gab_t=exp(-PI4 * ((t-(trans+1))/scale).^2 )';
         %             gab_f= exp(-pi*(scale/dimBase*(f-(freq+1))).^2);
         
         %          gab_t=exp(-PI4 * ((t-(trans+0.5))/scale).^2 )'; %popr. HAK
         %          gab_f= exp(-pi*(scale/dimBase*(f-(freq+0.5))).^2); %popr. HAK
         %          
         %          wig_gab=kron(gab_t, gab_f);
         %          NORM=sum(sum(wig_gab));
         %          if NORM~=0
         %             wig=wig+modulus.^2*wig_gab/NORM;    
         %          end
         
         
         %liczenie w osiach TxF
         %          g1=(BI_SQRT_PI/scale)*(tlow - trans);
         %          g2=(BI_SQRT_PI/scale)*(tupp - trans);
         %          hgab_t=(scale/4)*(erf(g2)-erf(g1))';
         %          
         %          gf1=(SQRT_PI*scale/dimBase)*(flow - freq);
         %          gf2=(SQRT_PI*scale/dimBase)*(fupp - freq);
         %          hgab_f=0.5*(dimBase/scale)*(erf(gf2)-erf(gf1));
         %          
         %          wig_gabTF=kron(hgab_t, hgab_f);
         %          
         %          NORM_TF=sum(sum(wig_gabTF));
         %          if NORM_TF~=0
         %             wigTF=wigTF+modulus.^2*wig_gabTF/NORM_TF;    
         %          end
         
         %liczenie w osiach XxY
         g1=(BI_SQRT_PI/scale)*(tlowx - trans);
         g2=(BI_SQRT_PI/scale)*(tuppx - trans);
         %hgab_tx=(scale/4)*(erf(g2)-erf(g1))';
         hgab_tx=(erf(g2)-erf(g1))'; %normowanie ponizej
         
         gf1=(SQRT_PI*scale/dimBase)*(flowy - freq);
         gf2=(SQRT_PI*scale/dimBase)*(fuppy - freq);
         %hgab_fy=0.5*(dimBase/scale)*(erf(gf2)-erf(gf1));
         hgab_fy=(erf(gf2)-erf(gf1));%normowanie ponizej
         
         wig_gabXY=kron(sparse(hgab_tx), sparse(hgab_fy));
         %           wig_gabXY=kron(hgab_tx, hgab_fy);
         %NORM_XY=sum(sum(wig_gabXY));
         
         g1_book=(BI_SQRT_PI/scale)*(book_min_t_points - trans);
         g2_book=(BI_SQRT_PI/scale)*(book_max_t_points - trans);
         gf1_book=(SQRT_PI*scale/dimBase)*(book_min_f_points - freq);
         gf2_book=(SQRT_PI*scale/dimBase)*(book_max_f_points - freq);

         %calka z czesci atomu ktora jest w granicach liczenia ksiazki MP
         NORM_XY = (erf(g2_book)-erf(g1_book))*(erf(gf2_book)-erf(gf1_book));
         if NORM_XY~=0
            wigXY=wigXY+modulus.^2*full(wig_gabXY)./NORM_XY;
         end
      end;
   else % dirac
      %       transTF=1+floor((trans-minT));% przesuwamy czestosc - floor -zeby pasowala do siatki wig, 
      % +1 bo w matlabie macierze maja indeksy od 1
      transXY=tx(1+floor((trans-minT)/DX));% przesuwamy czestosc - floor -zeby pasowala do siatki wig, 
      % +1 bo w matlabie macierze maja indeksy od 1
      
      if transXY>=1 & transXY<=Xsize
         %          wigTF(transTF,:)=wigTF(transTF,:)+(modulus^2)/Fsize; %normalizacja HAK
         wigXY(transXY,:)=wigXY(transXY,:)+(modulus^2)/Fsize; %normalizacja HAK
      end;
   end;
   if mod(k,hm) == 0, fprintf(1,'*');end
end;


%PJD
wigXY=wigXY';


fprintf(1,'\n');
xx=xx/Fsamp;
yy=yy*Fsamp/dimBase;







