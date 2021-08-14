function printFigure(name, f)
    if nargin < 2
        f = gcf;
    end
    print(f,'-depsc','-painters', [name, '.eps'])
    savefig(f, name);
end