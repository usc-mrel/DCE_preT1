function res=rmse(xref,x)
res=sqrt(abs(mean((x(:)-xref(:)).^2)));
end