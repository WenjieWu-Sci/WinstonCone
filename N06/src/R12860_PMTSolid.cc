#include "R12860_PMTSolid.hh"

#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include <cmath>

double R12860_PMTSolid::delta_torlerance = 1E-2*mm;

R12860_PMTSolid::R12860_PMTSolid(
        double R1, double R1z, double R2, double R3, double H, double H3)
    : m_R1(R1), m_R1z(R1z), m_R2(R2), m_R3(R3), m_H(H), m_H3(H3), numSide(4) {
    G4cout << __FILE__ << ":" << __LINE__ << std::endl;

    m_H_1_2 = m_H - m_R1z - m_H3;
    m_theta = atan((m_R2+m_R3)/(m_H_1_2));

    m_R1p = sqrt(pow(m_H_1_2,2)+pow(m_R2+m_R3,2)) - m_R2;

    G4cout << __FILE__ << ":" << __LINE__ << std::endl;
}

R12860_PMTSolid::~R12860_PMTSolid() {;}

G4VSolid*
R12860_PMTSolid::GetSolid(G4String solidname, double thickness) {
    // Calculate Parameter first
    double r1t = m_R1 + thickness;
    double r1zt = m_R1z + thickness;
    double r2t = m_R2 - thickness;
    double r3t = m_R3 + thickness;
    double r1zp = (m_R1p+thickness);
    double r4t = r1zp * sin(m_theta);
    double h1t = r1zp * cos(m_theta);
    double h2t = r2t * cos(m_theta);
    double h3t = m_H3 + thickness;

    // Show variables
    G4cout << "r1t: " << r1t/mm << " mm" << G4endl;
    G4cout << "r1zt: "<< r1zt/mm<< " mm" << G4endl;
    G4cout << "r2t: " << r2t/mm << " mm" << G4endl;
    G4cout << "r3t: " << r3t/mm << " mm" << G4endl;
    G4cout << "r1zp: "<< r1zp/mm<< " mm" << G4endl;
    G4cout << "r4t: " << r4t/mm << " mm" << G4endl;
    G4cout << "h1t: " << h1t/mm << " mm" << G4endl;
    G4cout << "h2t: " << h2t/mm << " mm" << G4endl;
    G4cout << "h3t: " << h3t/mm << " mm" << G4endl;
    // Construct the PMT Solid
    // * PART 1
    // G4Sphere* pmttube_solid_sphere = new G4Sphere(
    //                                         solidname+"_1_Sphere",
    //                                         0*mm, // R min
    //                                         r1t, // R max
    //                                         0*deg, // Start Phi
    //                                         360*deg, // Delta Phi
    //                                         0*deg, // Start Theta
    //                                         180*deg  // Delta Theta
    //                                         );
    G4Ellipsoid* pmttube_solid_sphere = new G4Ellipsoid(
                                            solidname+"_1_Ellipsoid",
                                            r1t, // pxSemiAxis
                                            r1t, // pySemiAxis
                                            r1zt // pzSemiAxis
                                            );
    // * PART 2
    G4Tubs* pmttube_solid_tube = new G4Tubs(
                                    solidname+"_2_Tube",
                                    0*mm,  /* inner */ 
                                    r4t+delta_torlerance, /* pmt_r */ 
                                    h2t/2+delta_torlerance, /* part 2 h */ 
                                    0*deg, 
                                    360*deg);
    G4Torus* pmttube_solid_torus = new G4Torus(
                                        solidname+"_2_Torus",
                                        0*mm,  // R min
                                        r2t+delta_torlerance, // R max
                                        (r2t+r3t), // Swept Radius
                                        -0.01*deg,
                                        360.01*deg);
    G4SubtractionSolid* pmttube_solid_part2 = new G4SubtractionSolid(
                                            solidname+"_part2",
                                            pmttube_solid_tube,
                                            pmttube_solid_torus,
                                            0,
                                            G4ThreeVector(0,0,-h2t/2)
                                            );
//    G4cout << "tube radius: " << (r4t+delta_torlerance)/mm << " mm" << G4endl;
//    G4cout << "tube height: " << (h2t+2*delta_torlerance)/mm << " mm" << G4endl;
//    G4cout << "torus radius: " << (r2t+delta_torlerance)/mm << " mm" << G4endl;
//    G4cout << "shift: " << (h2t/2)/mm << " mm" << G4endl;

    // * PART 3
    G4Tubs* pmttube_solid_end_tube = new G4Tubs(
                                    solidname+"_3_EndTube",
                                    0*mm,  /* inner */ 
                                    r3t+delta_torlerance, //21*cm/2, /* pmt_r */ 
                                    h3t/2+delta_torlerance, //30*cm/2, /* pmt_h */ 
                                    0*deg, 
                                    360*deg);

    // Join

    // * PART 1 + 2
    G4UnionSolid* pmttube_solid_1_2 = new G4UnionSolid(
                                            solidname+"_1_2",
                                            pmttube_solid_sphere,
                                            pmttube_solid_part2,
                                            0,
                                            G4ThreeVector(0, 0, -(h1t+h2t/2))
                                            );
    //return pmttube_solid_1_2;
    // * PART 1+2 + 3
    G4UnionSolid* pmttube_solid_1_2_3 = new G4UnionSolid(
                                            solidname,
                                            pmttube_solid_1_2,
                                            pmttube_solid_end_tube,
                                            0,
                                            G4ThreeVector(0,0, 
                                                -(m_H_1_2+h3t*0.50))
                                            );
    return pmttube_solid_1_2_3;
}

/*******************************************************************************
G4VSolid*
R12860_PMTSolid::GetCoverTop(G4String solidname, 
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double thickness) {

    G4VSolid* pmt_outer_solid = GetSolid(
                                    solidname+"_pmt_outer",
                                    thickness
                                    );

    G4double zPlane[2];
    G4double rInner[2] = {0, 0};
    G4double rOuter[2];

    zPlane[0] = c_h1;
    zPlane[1] = c_h2;

    rOuter[0] = c_r1;
    rOuter[1] = c_r2;

    G4Polyhedra* mHexagon = new G4Polyhedra(solidname+"_hexagon",
            0.*deg,         // statr phi
            360.*deg,       // total phi
            numSide,              // numSide   
            2,              // numZPlanes    
            zPlane,
            rInner,
            rOuter );

    G4VSolid* cover_solid_part_1 = new G4SubtractionSolid(solidname+"_part_1",
                                    mHexagon,
                                    pmt_outer_solid,
                                    0,
                                    G4ThreeVector());

    return cover_solid_part_1;


}

G4VSolid*
R12860_PMTSolid::GetCoverBottom(G4String solidname, 
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double thickness_out,
                            G4double thickness_in) {
    G4VSolid* pmt_outer_solid = GetSolid(
                                    solidname+"_pmt_bottom_outer",
                                    thickness_out
                                    );
    G4VSolid* pmt_inner_solid = GetSolid(
                                    solidname+"_pmt_bottom_inner",
                                    thickness_in
                                    );

    G4double zPlane[2];
    G4double rInner[2] = {0, 0};
    G4double rOuter[2];

    zPlane[0] = c_h2;
    zPlane[1] = c_h3;

    rOuter[0] = c_r2;
    rOuter[1] = c_r3;


    G4Polyhedra* mHexagon = new G4Polyhedra(solidname+"_hexagon_bottom",
            0.*deg,         // statr phi
            360.*deg,       // total phi
            numSide,        // numSide   
            2,              // numZPlanes    
            zPlane,
            rInner,
            rOuter );

    G4VSolid* cover_outer_all = new G4UnionSolid(solidname+"_bottom_outer_all",
                                                    pmt_outer_solid,
                                                    mHexagon,
                                                    0,
                                                    G4ThreeVector());

    G4VSolid* cover_outer_inner = new G4SubtractionSolid(solidname+"_bottom_outer_inner",
                                                        cover_outer_all,
                                                        pmt_inner_solid,
                                                        0,
                                                        G4ThreeVector());

    G4double zPlaneBigger[2];
    G4double rInnerBigger[2] = {0, 0};
    G4double rOuterBigger[2];

    zPlaneBigger[0] = c_h2;
    zPlaneBigger[1] = -(c_h2+m_H); // FIXME

    rOuterBigger[0] = c_r2*2;
    rOuterBigger[1] = c_r3*2;

    G4Polyhedra* mHexagonBigger = new G4Polyhedra(solidname+"_hexagon_bottom_bigger",
            0.*deg,         // statr phi
            360.*deg,       // total phi
            numSide,        // numSide   
            2,              // numZPlanes    
            zPlaneBigger,
            rInnerBigger,
            rOuterBigger);

    G4IntersectionSolid* cover_bottom = new G4IntersectionSolid(solidname,
                                                            cover_outer_inner,
                                                            mHexagonBigger,
                                                            0,
                                                            G4ThreeVector());

    return cover_bottom;
}

G4VSolid* 
R12860_PMTSolid::GetCoverSolid(G4String solidname, 
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double thickness) {
    // PART I
    G4double zPlane[3];
    G4double rInner[3] = {0, 0, 0};
    G4double rOuter[3];

    zPlane[0] = c_h1;
    zPlane[1] = c_h2;
    zPlane[2] = c_h3;

    rOuter[0] = c_r1;
    rOuter[1] = c_r2;
    rOuter[2] = c_r3;

    G4Polyhedra* mHexagon = new G4Polyhedra(solidname+"_hexagon_whole",
            0.*deg,         // statr phi
            360.*deg,       // total phi
            numSide,              // numSide   
            3,              // numZPlanes    
            zPlane,
            rInner,
            rOuter );

    // PART II
    G4VSolid* pmt_outer_solid = GetSolid(
                                    solidname+"_pmt_bottom_outer_bigger",
                                    thickness
                                    );

    // Union
    G4VSolid* pmt_cover_outer = new G4UnionSolid(
                                        solidname,
                                        pmt_outer_solid,
                                        mHexagon,
                                        0,
                                        G4ThreeVector());

    return pmt_cover_outer;

}


G4VSolid* 
R12860_PMTSolid::GetCoverTopCurve(G4String solidname,
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double center2equator,
                            G4double ball_r) {

    // PART I Sphere
    G4Sphere* top_solid_sphere = new G4Sphere(
                                            solidname+"_top_1_Sphere",
                                            0*mm, // R min
                                            ball_r, // R max
                                            0*deg, // Start Phi
                                            360*deg, // Delta Phi
                                            0*deg, // Start Theta
                                            180*deg  // Delta Theta
                                            );

    // PART II Polyhedra
    G4double zPlane[2];
    G4double rInner[2] = {0, 0};
    G4double rOuter[2];

    zPlane[0] = c_h1;
    zPlane[1] = c_h2;

    rOuter[0] = c_r1;
    rOuter[1] = c_r2;

    G4Polyhedra* mHexagon = new G4Polyhedra(solidname+"_top_curve_hexagon",
            0.*deg,         // statr phi
            360.*deg,       // total phi
            numSide,              // numSide   
            2,              // numZPlanes    
            zPlane,
            rInner,
            rOuter );

    // Subtraction
    G4VSolid* cover_solid = new G4SubtractionSolid(solidname,
                                    mHexagon,
                                    top_solid_sphere,
                                    0,
                                    G4ThreeVector(0, 0, center2equator));

    return cover_solid;

}

G4VSolid* 
R12860_PMTSolid::GetCoverSolidWithTopCurve(G4String solidname, 
                            G4double c_h0,
                            G4double c_r0,
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double center2equator,
                            G4double ball_r,
                            G4double thickness) {

    // PART I: main cover
    G4VSolid* main_cover = GetCoverSolid(
                                solidname+"_main_cover",
                                c_h1,
                                c_r1,
                                c_h2,
                                c_r2,
                                c_h3,
                                c_r3,
                                thickness
                            );
    // PART II: top curve cover
    G4VSolid* top_curve_cover = GetCoverTopCurve(
                                solidname+"_main_top_cover",
                                c_h0,
                                c_r0,
                                c_h1,
                                c_r1,
                                center2equator,
                                ball_r);

    // Union
    G4VSolid* whole_cover = new G4UnionSolid(
                                solidname,
                                main_cover,
                                top_curve_cover,
                                0,
                                G4ThreeVector());

    return whole_cover;
}
*******************************************************************************/
