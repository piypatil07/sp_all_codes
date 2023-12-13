# sp_all_codes
Exp No.1
TITLE: To generate and analyze different types of analog signals and discrete sequences using Matlab
t=0:0.1:2*pi; % from 0 to 2pi with steps of 0.1 radians

f_x=sin(t); % for Sin Function
subplot(4,1,1); 
plot (t,f_x,'k'); 
title('Subplot 1: sin(t) Continuous')

subplot(4,1,2); 
stem (t,f_x, 'k', 'filled');
title('Subplot 2: sin(t) Discrete')
 

f_y=cos(t); 
subplot(4,1,3)
plot (t,f_y, 'k'); 
title('Subplot 3: cos(t) Continuous')

subplot(4,1,4)
stem (t,f_y, 'k', 'filled'); 
title('Subplot 4: cos(t) Discrete')

t= -1:0.01:1;
unitstep = t>=0;
subplot(3,2,1); % To plot the graphs side by side
plot (t,unitstep, 'k'); % To plot t vs f_x
title('Subplot 1: unitstep Continuous')
subplot(3,2,2); % To plot the graphs side by side

stem (t,unitstep, 'k','filled'); % To plot t vs unitstep
title('Subplot 2: unitstep Discrete') % For assigning the titles to the graph.

impulse = t==0; % for Impulse Function.
subplot(3,2,3); % To plot the graphs side by side

plot (t,impulse, 'k'); % To plot t vs impulse
title('Subplot 3: impulse Continuous') % For assigning the titles to the graph.
subplot(3,2,4); % To plot the graphs side by side

stem (t,impulse, 'k','filled'); % To plot t vs impulse
title('Subplot 4: impulse Discrete') % For assigning the titles to the graph.

ramp = t.* unitstep;
subplot(3,2,5); % To plot the graphs side by side
plot (t,ramp, 'k'); % To plot t vs ramp
title('Subplot 5: ramp Continuous')
subplot(3,2,6); % To plot the graphs side by side

stem (t,ramp, 'k','filled'); % To plot t vs ramp
title('Subplot 6: ramp Discrete') % For assigning the titles to the graph



 EXP no.2
 TITLE: To perform basic signal manipulation like time-shifting, time reversal, time scaling using Matlab

n=0:0.5:9;
y=0:0.5:9;

m1=n-3;
subplot(4,2,1);
stem(m1,y);
title("Time Shifting")

m2=n+3;
subplot(4,2,2);
stem(m2,y);
title("Time Shifting")

m3=-n;
subplot(4,2,3);
stem (m3,y);
title("Time Reversal")

m4 = 3*n;
subplot(4,2,4);
stem(m4,y);
title("Time Scaling")

m5=(n-3)/3;
subplot(4,2,5);
stem(m5,y);
title("Shifting and Scaling")

m6=(n+3)/4;
subplot(4,2,6);
stem(m6,y);
title("Shifting and Scaling")

m7=n/4;
subplot(4,2,7);
stem(m7,y); 
title("Scaling")



Exp No.3
TITLE: To perform sampling, up sampling and down sampling using Matlab


t=-50:0.1:50;
fm=0.08;
f_x=cos(2*fm*t);
subplot(2,2,1)
plot(t,f_x)
title("Original Input Signal")

%under Sampling fs<2fm
n=-2:2;
fs=0.08;
x= cos(2*pi*fm*(n/fs));
subplot(2,2,2);
stem(n,x);
title("Under Sampling");

%critical Sampling fs=2fm
n=-4:4;
fs=0.16;
x= cos(2*pi*fm*(n/fs));
subplot(2,2,3);
stem(n,x);
title("Critical Sampling");

%over sampling fs>2fm
n=-8:8;
fs=1.5;
x= cos(2*pi*fm*(n/fs));
subplot(2,2,4);
stem(n,x);
title("Over Sampling")














x=[1,2,3,4,5,6,7,8,9,10,11,12];

a=upsample(x,3);

b=downsample(a,3);

c=downsample(b,3);

n=1:12;

subplot(2,2,1);
stem(n,x);
title("Original Input Signal")


%upsampling

n=1:36;

subplot(2,2,2);

stem(n,a);

title("Up Sampling");

%downsampling

n=1:4;

subplot(2,2,3);

stem(n,c);

title("Down Sampling");





Exp.No 4
TITLE: To perform linear convolution  of two sequences using Matlab
g=[2,4,6,9,1];
h=[-1,3,5,-6,3,2,5];
y= conv(g,h);
n1=length(g);
n2=length(h);
n=1:(n1+n2-1);
g=[g,zeros(1,n2-1)];
h=[h,zeros(1,n1-1)];
m=n1+n2-1;
z=zeros(1,m);
for i=1:m
 for j=1:i
 z(i)=z(i)+g(j)*h(i-j+1);
 end
end
disp(z)
subplot(2,1,1)
stem(1:m,z)
subplot(2,1,2)
stem(n,y)



exp no 5
TITLE: To perform auto-correlation using Matlab
N=1024;
f1=1;
fs=200;
n=0:N-1;
x=sin(2*pi*f1*n/fs);
t=(1:N)*(1/fs);
subplot(2,1,1);
plot(t,x);
xlabel('time')
ylabel('Aqmplitude')
grid;

Rxx=xcorr(x);
subplot(2,1,2);
plot(Rxx)
grid;

title ('auto corelation')
xlabel('lags')
ylabel('Auto corelation')






exp no.6
TITLE: To perform cross-correlation using Matlab
N=1024;
f1=1;
fs=200;
n=0:N-1;
x=sin(2*pi*f1*n/fs);
y=0.2*randn(1,N);
t=(1:N)*(1/fs);
subplot(3,1,1);
plot(t,x);
xlabel('Time');
ylabel('Amplitude');
grid;
subplot(3,1,2);
plot(t,y);
xlabel('Time')
ylabel('Amplitude')
grid;
Rxx=xcorr(x,y);
subplot(3,1,3);
plot(Rxx)
grid;
title ('Cross Correlation')
xlabel('lags')
ylabel('Cross Correlation')

