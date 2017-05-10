 clear all; 
rng(0);

nbits_per_channel = 4096;
nframes = 1000;
nerr = 100;

rsse=[4,1];
channel_epsilon = 0.3;
pr_channel = [1 2 1];
umat = [1 1; 1 -1];
lmat = [1/(1+channel_epsilon) 0; 0 1/(1-channel_epsilon)];
cmat = [1 channel_epsilon; channel_epsilon 1];
xin = [1 1 -1 -1;1 -1 -1 1];
zin = umat*xin;
trell = poly2bcjrtrellis(pr_channel, zin);
trellis=poly2rssetrellis(pr_channel, zin, rsse, trell);

% for s=1:3,
% x1=[0 0 0 0;0 0 0 0];
% [y1,zf]=filter(pr_channel,1,x1,[],2);
% end
snrvec = (9.5:0.5:12);
bervec = zeros(length(snrvec), 1);
bervec1 = zeros(length(snrvec), 1);

for i=1:length(snrvec),
    snr = snrvec(i);
    frame = 1;
    errc = 0;
    %e=zeros(nframes,1);
    for s=1:3,
        x1=[2 2 2 2;0 0 0 0];
        [y1,zf]=filter(pr_channel,1,x1,[],2);
    end
    while((frame<=nframes) && (errc < nerr))
        sigma2 = sum(pr_channel.^2)/(2*10.^((snr)/10));
        % transmitted symbols
        msg = randi(2, 2, nbits_per_channel)-1;
        %msg=rand(2, nbits_per_channel)>0.5;
        x = (2*msg - 1);
        z = umat * x;
     %npn=sqrt(sigma2).*randn(2, nbits_per_channel);
        npn = (lmat * umat * sqrt(sigma2)) *randn(2, nbits_per_channel);
        %rpn = lmat*umat*(cmat*filter(pr_channel, 1, x, zf, 2) + npn);
        rpn=filter(pr_channel, 1, z, zf, 2) + npn;
        y = vitdecoder(rpn, trell, channel_epsilon);
        y1=rssevitdecoder2(rpn, trellis, channel_epsilon,rsse);
        
        zout = zin(:, y);
        zout_rsse=zin(:,y1);
        xout =(umat)*zout;
        xout_rsse=(umat)*zout_rsse;
        msgout = (xout>0);
        rateerr = sum(sum(msg~=msgout))/(numel(msg));
        
         msg_rsse=(xout_rsse>0);
         rateerr1= sum(sum(msg~=msg_rsse))/(numel(msg));
         bervec1(i) = bervec1(i) + rateerr1;
        
%        if rateerr >0
%        e(frame,1)=frame;
%        end
        
        bervec(i) = bervec(i) + rateerr;
        frame = frame + 1;
        errc = errc + (rateerr1>0);
         
           
    end

    bervec(i) = bervec(i)/(frame-1);
    bervec1(i) = bervec1(i)/(frame-1);
end

figure
semilogy(snrvec, bervec,'r');
hold on
semilogy(snrvec, bervec1,'go--');
grid on;