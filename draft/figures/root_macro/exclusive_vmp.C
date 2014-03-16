void exclusive_vmp() {

    const unsigned line_width = 3;
    gStyle -> SetLineWidth(line_width);

    TCanvas * c = new TCanvas ("c1", "Exclusive Vector Meson Production", 1200, 1200);

    // incoming/outgoing electron
    TLine * e1 = new TLine(0.1, 0.9, 0.5, 0.7);
    TLine * e2 = new TLine(0.5, 0.7, 0.9, 0.9);
    e1 -> Draw();
    e2 -> Draw();

    // incoming/outgoing proton
    TLine * p1 = new TLine(0.1, 0.1, 0.5, 0.3);
    TLine * p2 = new TLine(0.5, 0.3, 0.9, 0.1);
    p1 -> Draw();
    p2 -> Draw();

    // photon
    TCurlyLine *gamma = new TCurlyLine(0.5, 0.7, 0.5, 0.5);
    gamma -> SetWavy();
    gamma -> Draw();

    // scalar object
    TLine * pom = new TLine(0.5, 0.3, 0.5, 0.5);
    pom -> SetLineStyle(2);
    pom -> Draw();

    // outgoing rho
    TLine * rho = new TLine (0.5, 0.5, 0.9, 0.5);
    rho -> Draw();

    // rho vertex
    TEllipse * vertex_rho = new TEllipse (0.5, 0.5, 0.015);
    vertex_rho -> SetFillColor(1);
    vertex_rho -> Draw();

    // proton vertex
    TEllipse * vertex_proton = new TEllipse (0.5, 0.3, 0.015);
    vertex_proton -> SetFillColor(1);
    vertex_proton -> Draw();

    // text labels
    TLatex text;

//     text.DrawLatex(0.1, 0.93,"e^{#pm}");
//     text.DrawLatex(0.9, 0.93,"e^{#pm}'");
    text.DrawLatex(0.1, 0.93,"h_{1}");
    text.DrawLatex(0.9, 0.93,"h_{1}'");


    text.DrawLatex(0.1, 0.05,"h_{2}");
    text.DrawLatex(0.9, 0.05,"h_{2}'");

    text.DrawLatex(0.9, 0.53,"V");

    c -> Print ("exclusive_vmp.eps");
}
