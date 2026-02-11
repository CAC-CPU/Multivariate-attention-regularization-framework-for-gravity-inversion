function main_inversion()
dataObs = load('data.txt');
mTrue = load('m_true.txt');
inversion = set_inversion_params(dataObs);
A = generate_forward_operator(inversion);
boundary = draw_model_boundary(inversion,mTrue);
switch inversion.kernelMode
    case 1
        inversion = trandition_reweight_conjugate_gradient_V1(inversion,A);
    case 2
        inversion = trandition_reweight_conjugate_gradient_V2(inversion,A);
    case 3
        inversion = teacher_reweight_conjugate_gradient(inversion,A);
    case 4
        inversion = modified_reweight_conjugate_gradient(inversion,A);
    case 5
        inversion = confidence_reweight_conjugate_gradient(inversion,A);
    case 6
        inversion = new_stabilizer_reweight_conjugate_gradient(inversion,A);
        plot_stablizers(inversion,boundary);
end
plot_inversion_params(inversion,boundary);
end