%R-state drive input: 
%to be used in conjunction with the
%Flip-flop switch simulation initializer script

function f = forceramp(t)
if t<100 ;
    f=0;
else
f = (t-100)*0.06;
end