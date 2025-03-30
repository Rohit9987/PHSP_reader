#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

/* code borrowed from EGSnrc henhouse IAEA record.

reads the phsp file.
need to modify the code to get the fluence at point defined
and then finally raytrace.
*/



#pragma pack(push, 1)  // Ensure no padding is added by the compiler

struct IAEA_Particle {
    short particle_type;       // Type of particle (Integer*2)
    float energy;              // Kinetic energy (MeV)
    float x, y, z;            // Position (cm)
    float u, v, w;            // Direction cosines
    float statistical_weight;  // Particle statistical weight
    char sign_of_w;            // Sign of W (Logical*1)
    char is_new_history;       // Is new history (Logical*1)
    int integer_extra[2];      // Extra integer storage
    float float_extra[2];      // Extra float storage
};

#pragma pack(pop)  // Restore default padding

int main() {
    std::ifstream file("../Varian_TrueBeam6MV_01.IAEAphsp", std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return 1;
    }

    IAEA_Particle particle;
    char particle_type_raw;

    std::cout << "Reading first 10 particles:\n";
    for (int i = 0; i < 100; i++) {
        // Read particle type (1 byte)
        if (!file.read(&particle_type_raw, sizeof(char))) break;
        
        // Extract sign of W from particle type
        short is = (particle_type_raw < 0) ? -1 : 1;
        particle.particle_type = std::abs(particle_type_raw);

        // Read energy (1 float)
        if (!file.read(reinterpret_cast<char*>(&particle.energy), sizeof(float))) break;
        
        // Determine if it's a new history
        particle.is_new_history = (particle.energy < 0) ? 1 : 0;
        particle.energy = std::fabs(particle.energy);
        
        // Read position (X, Y, Z) and direction cosines (U, V)
        if (!file.read(reinterpret_cast<char*>(&particle.x), sizeof(float) * 5)) break;

        
        // Compute W direction cosine
        double aux = (particle.u * particle.u + particle.v * particle.v);
        if (aux <= 1.0) {
            particle.w = static_cast<float>(is * sqrt(1.0 - aux));
        } else {
            aux = sqrt(aux);
            particle.u /= static_cast<float>(aux);
            particle.v /= static_cast<float>(aux);
            particle.w = 0.0f;
        }

        // Print particle information
        std::cout << "Particle " << i + 1 << ":\n";
        std::cout << "  Position: (" << particle.x << ", " << particle.y << ", " << particle.z << ") cm\n";
        std::cout << "  Direction: (" << particle.u << ", " << particle.v << ", " << particle.w << ")\n";
        std::cout << "  Energy: " << particle.energy << " MeV\n";
        std::cout << "  Type: " << particle.particle_type << "\n";
        std::cout << "  New history: " << (particle.is_new_history ? "Yes" : "No") << "\n\n";
    }

    file.close();
    return 0;
}

