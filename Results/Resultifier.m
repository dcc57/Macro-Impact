clear

format long
i = 1;
%for i = 1 : 17
    A = strcat('PreDisp',int2str(i));
    load(A,'PreDisp');
    A = PreDisp;
%end