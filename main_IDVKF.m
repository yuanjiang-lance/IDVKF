%
% A numerical example of IDVKF for dispersive signal decomposition
%
% Author: Yuan JIANG
% Time: 2026-06-06

clear, clc, close all;
load data_5dB;
%%
window = 128;
delta = 50;
beta1 = 1e-4;
bw = T/100;
tolout = 0.1;
r0 = [5e3, 5e2, 5e3];
tol = 1e-5;
maxit = 300;
betaC = 1e-7;
winLen = 1500;
%% extract ridges
[tindex,~] = ridgeDetectMult_F(SigFn, Fs, delta, beta1, bw, Nt, window, tolout);
figure(101), hold on, box on;
set(gca, 'LooseInset', [0,0,0,0]);
set(gcf,'Position',[50,340.2,400,350]);
set(gcf,'Color','w');
plot(t(tindex),f,'linewidth',1.5);
axis([0 T 0 Fs/2]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('GD before RPRG');
set(gca,'FontSize',13)
set(gca,'linewidth',1);
set(gca,'FontName','Times New Roman');

%% RPRG
thrt = length(t) / 40;
[tindexnew,~] = RPRG(tindex, thrt);
figure(201)
set(gca, 'LooseInset', [0,0,0,0]);
set(gcf,'Position',[450,340.2,400,350]);
set(gcf,'Color','w');
plot(t(tindexnew),f,'linewidth',1.5);
axis([0 T 0 Fs/2]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('GD after RPRG');
set(gca,'FontSize',13)
set(gca,'linewidth',1);
set(gca,'FontName','Times New Roman');

%%
iniGD = IFsmooth(t(tindexnew), 1e-5);
ID = zeros(3, 1);
[~, ID(1)] = min([norm(GD(1,:) - iniGD(1,:)), norm(GD(1,:) - iniGD(2,:)), norm(GD(1,:) - iniGD(3,:))]);
[~, ID(2)] = min([norm(GD(2,:) - iniGD(1,:)), norm(GD(2,:) - iniGD(2,:)), norm(GD(2,:) - iniGD(3,:))]);
[~, ID(3)] = min([norm(GD(3,:) - iniGD(1,:)), norm(GD(3,:) - iniGD(2,:)), norm(GD(3,:) - iniGD(3,:))]);
iniGD = iniGD(ID, :);
[SigFest, GDest, ~, r, SigFiter, GDiter] = IDVKF(SigFn, T, iniGD, r0, betaC, winLen, tol, maxit);
for k = 1: size(SigFest, 1)
    DFs(k, :) = [SigFest(k, :), conj(fliplr(SigFest(k, 2:ceil(Nt/2))))];
    SigTest(k, :) = real(ifft(DFs(k, :)));
end
figure(301), hold on, box on;
set(gca, 'LooseInset', [0,0,0,0]);
set(gcf,'Position',[1300,340.2,400,350]);
set(gcf,'Color','w');
plot(GDest,f,'b','linewidth',1.5);
plot(GD,f,'r--','linewidth',1.5)
axis([0 T 0 Fs/2]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Estimated GD');
set(gca,'FontSize',13)
set(gca,'linewidth',1);
set(gca,'FontName','Times New Roman');
%%
figure(401); clf
set(gcf,'Color','w');
set(gcf, 'Position', [120 700 270 270]);
set(gcf,'Color','w');
tlF = tiledlayout(3,1,'TileSpacing','none','Padding','none');

ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
plot(ax1, f, real(SigFest(1,:)), 'b-',  'LineWidth',0.5);
plot(ax1, f, real(SigFest(1,:) - SigF_mode(1,:)), 'r--','LineWidth',0.5);
plot(ax1, f, abs(SigFest(1,:)), 'k--','LineWidth',0.5);
ylim([-2.5, 2.5])
ylabel(ax1,'C1');
set(ax1,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);
ax1.XTickLabel = [];
ax1.XLabel.String = ''; 

ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
plot(ax2, f, real(SigFest(2,:)), 'b-',  'LineWidth',0.5);
plot(ax2, f, real(SigFest(2,:) - SigF_mode(2,:)), 'r--','LineWidth',0.5);
plot(ax2, f, abs(SigFest(2,:)), 'k--','LineWidth',0.5);
ylim([-2.5, 2.5])
ylabel(ax2,'C2');
set(ax2,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);
ax2.XTickLabel = [];
ax2.XLabel.String = '';

ax3 = nexttile; hold(ax3,'on'); box(ax3,'on');
plot(ax3, f, real(SigFest(3,:)), 'b-',  'LineWidth',0.5);
plot(ax3, f, real(SigFest(3,:) - SigF_mode(3,:)), 'r--','LineWidth',0.5);
plot(ax3, f, abs(SigFest(3,:)), 'k--','LineWidth',0.5);
ylim([-2.5, 2.5])
ylabel(ax3,'C3');
xlabel(ax3,'Frequency (Hz)');
set(ax3,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);

linkaxes([ax1 ax2 ax3],'x');
set([ax1 ax2 ax3],'LooseInset',[0 0 0 0]);

%%
figure(410); clf
set(gcf,'Position',[420 500 270 270]); 
set(gcf,'Color','w');
tlT = tiledlayout(3,1,'TileSpacing','none','Padding','none');

ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
plot(ax1, t, SigTest(1,:), 'b-', 'LineWidth',0.5);
plot(ax1, t, SigTest(1,:) - SigT_mode(1,:), 'r--', 'LineWidth',0.5);
ylim([-0.13, 0.13])
ylabel(ax1,'C1');
set(ax1,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);
ax1.XTickLabel = [];
ax1.XLabel.String = '';

ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
plot(ax2, t, SigTest(2,:), 'b-', 'LineWidth',0.5);
plot(ax2, t, SigTest(2,:) - SigT_mode(2,:), 'r--', 'LineWidth',0.5);
ylim([-0.08, 0.08])
ylabel(ax2,'C2');
set(ax2,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);
ax2.XTickLabel = [];
ax2.XLabel.String = '';

ax3 = nexttile; hold(ax3,'on'); box(ax3,'on');
plot(ax3, t, SigTest(3,:), 'b-', 'LineWidth',0.5);
plot(ax3, t, SigTest(3,:) - SigT_mode(3,:), 'r--', 'LineWidth',0.5);
ylim([-0.13, 0.13])
ylabel(ax3,'C3');
xlabel(ax3,'Time (s)');
set(ax3,'FontName','Times New Roman','FontSize',11,'LineWidth',0.6);

linkaxes([ax1 ax2 ax3],'x');
set([ax1 ax2 ax3],'LooseInset',[0 0 0 0]);
