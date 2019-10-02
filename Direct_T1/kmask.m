nx = 240; ny = 120;
c = sqrt(abs((nx/2)^2 - (ny/2)^2));

k_mask = ones(nx, ny);
for inx = 1:nx
    for iny = 1:ny
        if iny == 60
        end
        dist = sqrt(inx.^2 + (iny - c).^2) + sqrt(inx.^2 + (iny + c).^2);
        if dist > inx
            k_mask(inx, iny) = 0;
        end
    end
end

figure;
imagesc(k_mask)