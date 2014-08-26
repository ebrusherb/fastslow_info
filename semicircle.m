function semicircle(center,radius,color,lw,semi)
hold on
if semi
    radians=-3*pi/4:.01:pi/4;
    xcirc=radius*cos(radians)+center(1);
    ycirc=radius*sin(radians)+center(2);
    fill(xcirc,ycirc,color,'EdgeColor','none','LineWidth',lw)
    radians=pi/4:.01:5*pi/4;
    xcirc=radius*cos(radians)+center(1);
    ycirc=radius*sin(radians)+center(2);
    fill(xcirc,ycirc,'k','EdgeColor','none','LineWidth',lw)
else 
    radians=-3*pi/4:.01:5*pi/4;
    xcirc=radius*cos(radians)+center(1);
    ycirc=radius*sin(radians)+center(2);
    fill(xcirc,ycirc,color,'EdgeColor','none','LineWidth',lw)
end
radians=-3*pi/4:.01:5*pi/4;
xcirc=radius*cos(radians)+center(1);
ycirc=radius*sin(radians)+center(2);
plot(xcirc,ycirc,'k','LineWidth',lw)

