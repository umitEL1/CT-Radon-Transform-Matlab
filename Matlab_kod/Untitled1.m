clc
close all;
clear all;
disp("Boyutu 256*256, kaynak derecesi 0-180 olan phantom modeli i�in ���n say�s�n� giriniz. ");
disp("I��n say�s� netlik oran�n� do�rudan etkilemektedir. ");
t=input('I�IN SAYISI = ');                          % ���n say�s�      

resim = phantom(256);                               % resim tan�m
[SA,ST] = size(resim);                              % resim boyutu
sa = round(SA/2);                                   % sat�r yar� boyut
st = round(ST/2);                                   % s�tun yar� boyut
kosegen = ceil(sqrt(power(SA,2)+power(ST,2)));      % k��egen uzunlu�u(kaynak uzunlu�u)
k = ceil(kosegen/2);                                % k��egen yar�s�
kaynak_araligi = (kosegen/(t-1));   % t ���n var (t-1) bo�luk,iki k��egen aras� mesafe

izdusum = zeros(t,180);             % s�f�rlardan olu�an izd���m matrisi
resim(:,ST+1:ST+4)=0;      %remin s�n�rkar�na s�f�r ekledim
resim(SA+1:SA+4,:)=0;      %baz� tan�ms�z olma durumlar�n� engelledim

subplot(2,2,1),imshow(resim);  %orjinal resmi �izdirdi
title('ORJ�NAL RES�M');

for q = 1:45  
     a = tand(90+q);                              % e�im
     v=t;
     for r = -k:kaynak_araligi:k                  % ���n noktalar�
         v=v-1;                                   % indis belirleme
         b = r / sind(q);                         % x=0 iken y eksenini kesti�i nokta
         bitis = round(min((-a*st+b),sa));        % ��kt��� nokta y yi baz alarak
         baslangic = round(max((a*st+b),-sa));    % girdi y baz al�narak
         if (bitis > st || bitis < -st || baslangic > st || baslangic < -st)
             bitis = 0 ;                   % ���n�n resime �arpmad��� durumlar
             baslangic = 0 ;
         end
         if(baslangic < bitis)             % resime �arpmayan ���nlar� i�lemden engelleme
             for y = baslangic :(bitis-1)  % y de�erini ba�lang��tan biti�e sayd�rd�m
                 x=(y-b)/a;                % y=ax+b x notas�n� saptad�m
                 R=1/cosd(q); %y'nin 1birim artmas� durumumda ���n�n x deki toplam a��rl��� 
                 x_tam= floor(x); % x i floor komutu ile bir alt tam say�ya yuvarlad�m
                 parca_x= x-x_tam; % bu iki degeri birbirinden ��kararak k�d�rat� buldum
                 r1=parca_x/cosd(q); % kusurattan faydalanarak x deki ilk a��rl��� buldum
                 r2=R-r1;  % toplam a��rl�ktan r1 i ��kararak 2. a��rl��� hesaplad�m
                 if r1>R    
                     r1=R;  % r1 in x ekseninde tek bir piksen ge�ti�i �zel durumlar 
                     r2=0;
                 end
                 while r1 > 1.5 || r2>1.5         % algoritmam�n dogu �al���p �al��mad���n� 
                     disp("hataaaaaaaaaaaaaaaaa");%test ettim r1 ve r2 de�erleri
                 end                              %0-k�k2 aras� de�erler olml�d�r
                 izdusum(v,q)= izdusum(v,q)+r1*resim(y+sa+1 , x_tam+st+2);%a��rl�klar� pikseller ile �arpt�m
                 izdusum(v,q)= izdusum(v,q)+r2*resim(y+sa+1 , x_tam+st+1);% izd�s�m matrisinde uygun yere arad�m 
             end
         end 
     end
     subplot(2,2,2),imagesc(izdusum);       %1-45 aras� sinogram� �izdirdim
     xlabel("A�I = " + q); ylabel("TOPLAM I�IN SAYINIZ = "+t);%eksen adland�rma
     title('�ZD���M MATR�S� = A�IRLIK * P�KSEL'); % ba�l�k
     drawnow %anl�k �izim �zelli�i
     colormap(gray); % g�ri seviye �izim
end     
for q = 46:90  
     a = tand(90+q);
     v=t;
     for r = -k:kaynak_araligi:k
         v=v-1;
         b = r / sind(q);
         bitis = round(min((-st-b)/a,st));   
         baslangic = round(max((st-b)/a,-st));
         if (bitis>st || bitis<-st || baslangic>st || baslangic<-st)
             bitis = 0 ;   
             baslangic = 0 ;
         end
         if(baslangic<bitis)
             for x=baslangic:(bitis-1)
                 y=a*x+b;
                 R=1/sind(q);
                 y_tam= ceil(y);
                 parca_y= y_tam-y;
                 r1=parca_y/sind(q);
                 r2=R-r1; 
                 if r1>R
                     r1=R;
                     r2=0;
                 end
                 while r1 > 1.5 || r2>1.5 || r1<0 || r2<0
                     disp("hataaaaaaaaaaaaaaaaa");
                 end
                 izdusum(v,q)= izdusum(v,q)+r1*resim(y_tam+sa+1, x+st+1);
                 izdusum(v,q)= izdusum(v,q)+r2*resim(y_tam+sa+2  , x+st+1);   
             end
         end
     end
     subplot(2,2,2),imagesc(izdusum);
     xlabel("A�I = " + q); ylabel("TOPLAM I�IN SAYINIZ = "+t);
     title('�ZD���M MATR�S� = A�IRLIK * P�KSEL');
     drawnow
     colormap(gray);
end   
for q = 91:135  
     a = tand(q-90);
     v=t;
     for r = -k:kaynak_araligi:k
         v=v-1;
         b = r / sind(q);
         bitis = round(min((st-b)/a,st));   %y=a*x+b girdi
         baslangic = round(max((-st-b)/a,-st));
         if (bitis>st ||bitis<-st ||baslangic>st || baslangic<-st)
             bitis= 0 ;   
             baslangic = 0 ;
         end
         if(baslangic<bitis)
             for x=baslangic:(bitis-1)
                 y=a*x+b;
                 R=1/sind(q);
                 y_tam= ceil(y);
                 parca_y= y_tam-y;
                 r1=parca_y/sind(q);
                 r2=R-r1; 
                 if r1>R
                     r1=R;
                     r2=0;
                 end
                 while r1 > 1.5 || r2>1.5 || r1<0 || r2<0
                     disp("hataaaaaaaaaaaaaaaaa");
                 end
                 izdusum(v,q)= izdusum(v,q)+r1*resim(y_tam+sa+1, x+st+1);
                 izdusum(v,q)= izdusum(v,q)+r2*resim(y_tam+sa+2  , x+st+1);    
             end
         end   
     end
     subplot(2,2,2),imagesc(izdusum);
     xlabel("A�I = " + q); ylabel("TOPLAM I�IN SAYINIZ = "+t);
     title('�ZD���M MATR�S� = A�IRLIK * P�KSEL');
     drawnow
     colormap(gray);
end 
for q = 136:180  
     a = tand(90+q);
     v=t;
     for r = -k:kaynak_araligi:k
         v=v-1;
         b = r / sind(q);
         bitis = round(min((a*st+b),sa));   %y=a*x+b
         baslangic = round(max((-a*st+b),-sa));
         if (bitis>st || bitis<-st ||baslangic>st ||baslangic<-st)
             bitis = 0 ;   
             baslangic = 0 ;
         end
         if(baslangic<bitis)
             for y=baslangic:(bitis-1)
                 x=(y-b)/a;
                 R=(1/cosd(q));
                 x_tam= ceil(x);
                 parca_x= x_tam - x;
                 r1=(parca_x/cosd(q));
                 r2=R-r1; 
                 if r1>R
                     r1=R;
                     r2=0;
                 end
                 while r1 > 1.5 || r2>1.5
                     disp("hataaaaaaaaaaaaaaaaa");
                 end
                 izdusum(v,q)= izdusum(v,q)+(-r1*resim(y+sa+1, x_tam+st+2));
                 izdusum(v,q)= izdusum(v,q)+(-r2*resim(y+sa+1, x_tam+st+1));   
             end
         end 
     end
     subplot(2,2,2),imagesc(izdusum);
     xlabel("A�I = " + q); ylabel("TOPLAM I�IN SAYINIZ = "+t);
     title('�ZD���M MATR�S� = A�IRLIK * P�KSEL');
     drawnow
     colormap(gray);
end

geriizdusum = zeros(t,t);          %s�f�rlardan olu�an bir geri izd�s�m matrisi tan�mlad�k
mid = floor(t/2)+1;
[x,y] = meshgrid(ceil(-t/2):ceil(t/2-1));
for i = 1:180
    rot = round(mid+ x*sind(271-i) +y*cosd(271-i));
    indis   = find((rot > 0) & (rot <= t));
    Coordinat = rot(indis);  
    geriizdusum(indis) = geriizdusum(indis) + izdusum(Coordinat,i)./180;
    subplot(2,2,3),imagesc( geriizdusum);                          %geri izd�s�m� �izdirdim
    xlabel("A�I = " + i); ylabel("TOPLAM I�IN SAYINIZ = "+t);      %eksenlere isim verdim
    title('GER� �ZD���M MATR�S� ');                                %ba�l�k
    drawnow                                                        %anl�k �izimler
    colormap(gray);                                                %gri ton
end

SG = fftshift(fft(geriizdusum')',2); %FFD geri izd���m�n fouriyerini ald�m
ramp = ones(t,1)*abs(-ceil((t-1)/2):floor((t-1)/2)); % rampa fitresi
FGI = ifft(ifftshift(SG.*ramp,2)')'; %filtrelenmi� geri izd�s�m =ters ffd(rampa*ffd)
subplot(2,2,4),imagesc(FGI);
title('F�LTRELENM�� GER� �ZD���M MATR�S�');
colormap(gray);
