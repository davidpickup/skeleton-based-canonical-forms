% Load the sample mesh.
load mesh.mat

% Compute a canonical form using our method.
canonical = skeletonCanonicalForm(mesh);
save canonical.mat canonical

% Display the results.
figure;
subplot(1,2,1);
visualiseMesh(mesh);
title('Original mesh');
subplot(1,2,2);
visualiseMesh(canonical);
title('Canonical form');
