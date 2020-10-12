clc
close all;
clear all;
disp("Boyutu 256*256, kaynak derecesi 0-180 olan phantom modeli için ýþýn sayýsýný giriniz. ");
disp("Iþýn sayýsý netlik oranýný doðrudan etkilemektedir. ");
t=input('IÞIN SAYISI = ');                          % ýþýn sayýsý      

resim = phantom(256);                               % resim taným
[SA,ST] = size(resim);                              % resim boyutu
sa = round(SA/2);                                   % satýr yarý boyut
st = round(ST/2);                                   % sütun yarý boyut
kosegen = ceil(sqrt(power(SA,2)+power(ST,2)));      % köþegen uzunluðu(kaynak uzunluðu)
k = ceil(kosegen/2);                                % köþegen yarýsý
kaynak_araligi = (kosegen/(t-1));   % t ýþýn var (t-1) boþluk,iki köþegen arasý mesafe

izdusum = zeros(t,180);             % sýfýrlardan oluþan izdüþüm matrisi
resim(:,ST+1:ST+4)=0;      %remin sýnýrkarýna sýfýr ekledim
resim(SA+1:SA+4,:)=0;      %bazý tanýmsýz olma durumlarýný engelledim

subplot(2,2,1),imshow(resim);  %orjinal resmi çizdirdi
title('ORJÝNAL RESÝM');

for q = 1:45  
     a = tand(90+q);                              % eðim
     v=t;
     for r = -k:kaynak_araligi:k                  % ýþýn noktalarý
         v=v-1;                                   % indis belirleme
         b = r / sind(q);                         % x=0 iken y eksenini kestiði nokta
         bitis = round(min((-a*st+b),sa));        % çýktýðý nokta y yi baz alarak
         baslangic = round(max((a*st+b),-sa));    % girdi y baz alýnarak
         if (bitis > st || bitis < -st || baslangic > st || baslangic < -st)
             bitis = 0 ;                   % ýþýnýn resime çarpmadýðý durumlar
             baslangic = 0 ;
         end
         if(baslangic < bitis)             % resime çarpmayan ýþýnlarý iþlemden engelleme
             for y = baslangic :(bitis-1)  % y deðerini baþlangýçtan bitiþe saydýrdým
                 x=(y-b)/a;                % y=ax+b x notasýný saptadým
                 R=1/cosd(q); %y'nin 1birim artmasý durumumda ýþýnýn x deki toplam aðýrlýðý 
                 x_tam= floor(x); % x i floor komutu ile bir alt tam sayýya yuvarladým
                 parca_x= x-x_tam; % bu iki degeri birbirinden çýkararak küdüratý buldum
                 r1=parca_x/cosd(q); % kusurattan faydalanarak x deki ilk aðýrlýðý buldum
                 r2=R-r1;  % toplam aðýrlýktan r1 i çýkararak 2. aðýrlýðý hesapladým
                 if r1>R    
                     r1=R;  % r1 in x ekseninde tek bir piksen geçtiði özel durumlar 
                     r2=0;
                 end
                 while r1 > 1.5 || r2>1.5         % algoritmamýn dogu çalýþýp çalýþmadýðýný 
                     disp("hataaaaaaaaaaaaaaaaa");%test ettim r1 ve r2 deðerleri
                 end                              %0-kök2 arasý deðerler olmlýdýr
                 izdusum(v,q)= izdusum(v,q)+r1*resim(y+sa+1 , x_tam+st+2);%aðýrlýklarý pikseller ile çarptým
                 izdusum(v,q)= izdusum(v,q)+r2*resim(y+sa+1 , x_tam+st+1);% izdüsüm matrisinde uygun yere aradým 
             end
         end 
     end
     subplot(2,2,2),imagesc(izdusum);       %1-45 arasý sinogramý çizdirdim
     xlabel("AÇI = " + q); ylabel("TOPLAM IÞIN SAYINIZ = "+t);%eksen adlandýrma
     title('ÝZDÜÞÜM MATRÝSÝ = AÐIRLIK * PÝKSEL'); % baþlýk
     drawnow %anlýk çizim özelliði
     colormap(gray); % gýri seviye çizim
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
     xlabel("AÇI = " + q); ylabel("TOPLAM IÞIN SAYINIZ = "+t);
     title('ÝZDÜÞÜM MATRÝSÝ = AÐIRLIK * PÝKSEL');
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
     xlabel("AÇI = " + q); ylabel("TOPLAM IÞIN SAYINIZ = "+t);
     title('ÝZDÜÞÜM MATRÝSÝ = AÐIRLIK * PÝKSEL');
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
     xlabel("AÇI = " + q); ylabel("TOPLAM IÞIN SAYINIZ = "+t);
     title('ÝZDÜÞÜM MATRÝSÝ = AÐIRLIK * PÝKSEL');
     drawnow
     colormap(gray);
end

geriizdusum = zeros(t,t);          %sýfýrlardan oluþan bir geri izdüsüm matrisi tanýmladýk
mid = floor(t/2)+1;
[x,y] = meshgrid(ceil(-t/2):ceil(t/2-1));
for i = 1:180
    rot = round(mid+ x*sind(271-i) +y*cosd(271-i));
    indis   = find((rot > 0) & (rot <= t));
    Coordinat = rot(indis);  
    geriizdusum(indis) = geriizdusum(indis) + izdusum(Coordinat,i)./180;
    subplot(2,2,3),imagesc( geriizdusum);                          %geri izdüsümü çizdirdim
    xlabel("AÇI = " + i); ylabel("TOPLAM IÞIN SAYINIZ = "+t);      %eksenlere isim verdim
    title('GERÝ ÝZDÜÞÜM MATRÝSÝ ');                                %baþlýk
    drawnow                                                        %anlýk çizimler
    colormap(gray);                                                %gri ton
end

SG = fftshift(fft(geriizdusum')',2); %FFD geri izdüþümün fouriyerini aldým
ramp = ones(t,1)*abs(-ceil((t-1)/2):floor((t-1)/2)); % rampa fitresi
FGI = ifft(ifftshift(SG.*ramp,2)')'; %filtrelenmiþ geri izdüsüm =ters ffd(rampa*ffd)
subplot(2,2,4),imagesc(FGI);
title('FÝLTRELENMÝÞ GERÝ ÝZDÜÞÜM MATRÝSÝ');
colormap(gray);
