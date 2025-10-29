// BallOnTriangles.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <vector>
#include <array>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) :x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return { x + o.x,y + o.y,z + o.z }; }
    Vec3 operator-(const Vec3& o) const { return { x - o.x,y - o.y,z - o.z }; }
    Vec3 operator*(double s) const { return { x * s,y * s,z * s }; }
    Vec3 operator/(double s) const { return { x / s,y / s,z / s }; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

static inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
static inline double len(const Vec3& a) { return sqrt(dot(a, a)); }
static inline Vec3 normalized(const Vec3& a) { double L = len(a); return L > 0 ? a / L : a; }

struct Tri {
    // indices of vertices (counter-clockwise)
    int a, b, c;
    // neighbors across edges (a-b), (b-c), (c-a); -1 means boundary
    int n_ab, n_bc, n_ca;
};

struct Mesh {
    vector<Vec3> V;
    vector<Tri>  T;
};

struct Bary {
    double u, v, w; // with respect to (a,b,c)
};

// Compute barycentric coordinates of p w.r.t. triangle (A,B,C)
static Bary barycentric(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C) {
    Vec3 v0 = B - A;
    Vec3 v1 = C - A;
    Vec3 v2 = P - A;
    double d00 = dot(v0, v0);
    double d01 = dot(v0, v1);
    double d11 = dot(v1, v1);
    double d20 = dot(v2, v0);
    double d21 = dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;
    return { u,v,w };
}

static Vec3 triNormal(const Vec3& A, const Vec3& B, const Vec3& C) {
    Vec3 n = normalized(cross(B - A, C - A));
    if (n.z < 0.0) n = n * -1.0; // alltid pek oppover
    return n;
}

// Project vector g to plane with normal n (remove normal component)
static Vec3 projectToPlane(const Vec3& g, const Vec3& n) {
    return g - n * dot(g, n);
}

// Clamp point P onto edge segment E0-E1
static Vec3 closestPointOnSegment(const Vec3& P, const Vec3& E0, const Vec3& E1) {
    Vec3 e = E1 - E0;
    double t = dot(P - E0, e) / max(1e-12, dot(e, e));
    t = max(0.0, min(1.0, t));
    return E0 + e * t;
}

int main() {
    ios::sync_with_stdio(false);

    // --- Manuell triangulering ---
    // Fem hjørner: fire i et kvadrat + ett senter. Vi lager 4 trekanter rundt senteret.
    Mesh M;
    M.V = {
	    {-2.0, -1.0, 0.4},   // v0 - venstre bunn
	    { 0.0, -1.0, 0.2},   // v1 - midt bunn (lavere)
	    {-2.0,  1.0, 0.4},   // v2 - venstre topp
	    { 0.0,  1.0, 0.0},   // v3 - midt topp
	    { 2.0, -1.0, 0.4},   // v4 - høyre bunn
	    { 2.0,  1.0, 0.4}    // v5 - høyre topp
    };

    // Indekser (mot klokka sett ovenfra)
    // T0, T1 = venstre plan
    // T2, T3 = høyre plan
    M.T = {
	    {0, 1, 2,  -1, 1, -1},   // T0
	    {1, 2, 3,   0, -1, 2},   // T1 (snudd rekkefølge)
	    {1, 3, 4,   1, 3, -1},   // T2 (snudd rekkefølge)
	    {4, 5, 3,   2, -1, -1}   // T3
    };

    for (size_t i = 0; i < M.T.size(); ++i) {
        auto t = M.T[i];
        Vec3 n = triNormal(M.V[t.a], M.V[t.b], M.V[t.c]);
        cout << "T" << i << " normal = (" << n.x << ", " << n.y << ", " << n.z << ")\n";
    }

    // Print triangulering og naboer
    cout << "Triangulering (4 trekanter) med naboer:\n";
    for (size_t i = 0; i < M.T.size(); ++i) {
        const auto& t = M.T[i];
        cout << "T" << i << ": (" << t.a << "," << t.b << "," << t.c << ")  "
            << "n_ab=" << t.n_ab << "  n_bc=" << t.n_bc << "  n_ca=" << t.n_ca << "\n";
    }
    cout << "\n";

    // --- Simulasjon ---
    // Tyngdeakselerasjon
    const Vec3 g = { 0.0, 0.0, -9.81 };

    // Friksjonskoeffisient (kinetisk, enkel modell – sett til 0 for ren rulling uten glidning)
    const double mu = 0.0;

    // Startposisjon: inne i T0 (bruk barysentriske for en trygg posisjon)
    int tri = 0;
    {
        const auto& t = M.T[tri];
        Vec3 A = M.V[t.a], B = M.V[t.b], C = M.V[t.c];
        double u = 0.2, v = 0.2, w = 0.6; // u+v+w=1
        Vec3 P = A * u + B * v + C * w;
        // Vi bruker P videre:
        // men vi legger alt i variabler utenfor blokka
    }

    // Sim-state
    Vec3 pos{-1.50, -0.5, 0.8};
	/*{
        const auto& t = M.T[tri];
        Vec3 A = M.V[t.a], B = M.V[t.b], C = M.V[t.c];
        double u = 0.2, v = 0.2, w = 0.6;
        pos = A * u + B * v + C * w;
    }*/
    Vec3 vel = { 0.0, 0.0, 0.0 };

    // Integrator
    const double dt = 0.01;
    const double T_end = 5.0;
    const double print_every = 0.1;
    double next_print = 0.0;

    auto triVerts = [&](int tIdx) { return array<Vec3, 3>{ M.V[M.T[tIdx].a], M.V[M.T[tIdx].b], M.V[M.T[tIdx].c] }; };
    auto triNeighbors = [&](int tIdx) { return array<int, 3>{ M.T[tIdx].n_ab, M.T[tIdx].n_bc, M.T[tIdx].n_ca }; };

    cout << fixed << setprecision(4);
    cout << "t(s)\ttri\tpos.x\tpos.y\tpos.z\t|v|\n";

    for (double t = 0.0; t <= T_end + 1e-9; t += dt) {

        auto Vtx = triVerts(tri);
        Vec3 A = Vtx[0], B = Vtx[1], C = Vtx[2];
        Vec3 n = triNormal(A, B, C);

        // Akselerasjon = gravitasjon projisert til planet - enkel friksjon mot bevegelsesretning
        Vec3 a_tangent = projectToPlane(g, n);
        Vec3 v_tangent = projectToPlane(vel, n);
        if (mu > 0.0 && len(v_tangent) > 1e-9) {
            Vec3 vhat = normalized(v_tangent);
            // Coulomb friksjon ~ mu * N, der N = |g·n| (per masseenhet)
            Vec3 friction = vhat * (-mu * fabs(dot(g, n)));
            a_tangent += friction;
        }

        // Semi-implisitt Euler
        vel += a_tangent * dt;
        pos += vel * dt;

        // Hold oss på aktuell trekant / flytt over kanter ved behov
        // Sjekk barysentriske etter oppdatert pos
        Bary bc = barycentric(pos, A, B, C);
        const double eps = -1e-8;

        if (bc.u >= eps && bc.v >= eps && bc.w >= eps) {
            // fortsatt inni – prosjektér pos på planet for numerisk stabilitet
            // pos_on_plane = A + (I - n n^T) * (pos - A)
            Vec3 AP = pos - A;
            pos = A + projectToPlane(AP, n);
        }
        else {
            // Finn mest negativ koeff – den kanten krysset vi
            double vals[3] = { bc.u, bc.v, bc.w };
            int mostNeg = 0;
            for (int i = 1; i < 3; ++i) if (vals[i] < vals[mostNeg]) mostNeg = i;

            // Edge-indeksene matcher nabo-rekkefølgen: 0->(a-b), 1->(b-c), 2->(c-a)
            int edge = -1;
            if (mostNeg == 0) edge = 2; // u<0 betyr utenfor ved kanten (b-c)? Nei, standard:
            // For barycentrics (u,v,w) with vertices (A,B,C):
            // u<0 => outside edge opposite A => edge (B,C) which corresponds to neighbor n_bc (index 1)
            // v<0 => outside (C,A) => neighbor n_ca (index 2)
            // w<0 => outside (A,B) => neighbor n_ab (index 0)
            if (mostNeg == 0) edge = 1; // u<0 => (B,C)
            if (mostNeg == 1) edge = 2; // v<0 => (C,A)
            if (mostNeg == 2) edge = 0; // w<0 => (A,B)

            auto NBs = triNeighbors(tri);
            int nb = NBs[edge];

            // Kanteendepunkter for denne kanten:
            Vec3 E0, E1;
            if (edge == 0) { E0 = A; E1 = B; }
            else if (edge == 1) { E0 = B; E1 = C; }
            else { E0 = C; E1 = A; }

            if (nb >= 0) {
                // Bytt til nabo-triangel
                tri = nb;
                auto Vtx2 = triVerts(tri);
                Vec3 A2 = Vtx2[0], B2 = Vtx2[1], C2 = Vtx2[2];
                Vec3 n2 = triNormal(A2, B2, C2);

                // Projiser posisjon til nytt plan
                double dist = dot(pos - A2, n2);
                pos -= n2 * dist;
                pos += n2 * 1e-4; // løfter den litt

                // Sjekk barysentriske og klem inn hvis nødvendig
                Bary bc2 = barycentric(pos, A2, B2, C2);
                if (bc2.u < 0 || bc2.v < 0 || bc2.w < 0) {
                    // klem til nærmeste punkt på triangelet
                    // (en enkel variant som sjekker kanter)
                    Vec3 Pab = closestPointOnSegment(pos, A2, B2);
                    Vec3 Pbc = closestPointOnSegment(pos, B2, C2);
                    Vec3 Pca = closestPointOnSegment(pos, C2, A2);

                    double dab = len(pos - Pab);
                    double dbc = len(pos - Pbc);
                    double dca = len(pos - Pca);
                    double dmin = min({ dab, dbc, dca });

                    if (dmin == dab) pos = Pab;
                    else if (dmin == dbc) pos = Pbc;
                    else pos = Pca;
                }

                // Oppdater hastighet til nytt plan
                vel = projectToPlane(vel, n2);
            }
            else {
                // Ytterrandskant: klamp pos til kanten og reflekter hastighet i planet rundt kanten (med litt demping)
                pos = closestPointOnSegment(pos, E0, E1);

                // Finn triangelplanet og kanttangent
                Vec3 nplane = n;
                Vec3 edgeDir = normalized(E1 - E0);
                // Decomposer v_tangent i langs kant og normal til kant (i planet)
                Vec3 vtan = projectToPlane(vel, nplane);
                Vec3 v_along = edgeDir * dot(vtan, edgeDir);
                Vec3 v_perp = vtan - v_along;

                double restitution = 0.2; // demp
                vel = v_along - v_perp * restitution; // "sprett" bort fra kanten
                // liten ekstra demping
                vel = vel * 0.98;
            }
        }

        if (t + 1e-12 >= next_print) {
            cout << ",\n("
                << pos.x << ", " << pos.y << ", " << pos.z << ")";
            next_print += print_every;
        }
    }

    return 0;
}
