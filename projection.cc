#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#pragma pack(push, 1)
struct IAEA_Particle {
    short particle_type;       // Type of particle
    float energy;              // Kinetic energy (MeV)
    float x, y, z;             // Position (cm)
    float u, v, w;             // Direction cosines
    float statistical_weight;  // Particle statistical weight
    char sign_of_w;            // Sign of W
    char is_new_history;       // Is new history
};
#pragma pack(pop)

int main() {
    std::ifstream file("../Varian_TrueBeam6MV_01.IAEAphsp", std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return 1;
    }

    IAEA_Particle particle;
    char particle_type_raw;
    int count_in_field = 0;
    int total_particles = 1000000; // Read 1 million particles
    float field_half_size = 5.0;   // Half of 10x10 cm^2 field
    float projection_plane = 100.0;

    for (int i = 0; i < total_particles; i++) {
        if (!file.read(&particle_type_raw, sizeof(char))) break;
        short is = (particle_type_raw < 0) ? -1 : 1;
        particle.particle_type = std::abs(particle_type_raw);

        if (!file.read(reinterpret_cast<char*>(&particle.energy), sizeof(float))) break;
        particle.is_new_history = (particle.energy < 0) ? 1 : 0;
        particle.energy = std::fabs(particle.energy);

        if (!file.read(reinterpret_cast<char*>(&particle.x), sizeof(float) * 5)) break;

        double aux = (particle.u * particle.u + particle.v * particle.v);
        particle.w = (aux <= 1.0) ? static_cast<float>(is * sqrt(1.0 - aux)) : 0.0f;

        // Project to z = 100 cm
        if (particle.w != 0) {
            float t = (projection_plane - particle.z) / particle.w;
            float x_proj = particle.x + particle.u * t;
            float y_proj = particle.y + particle.v * t;

            // Check if within 10x10 cm field
            if (std::abs(x_proj) <= field_half_size && std::abs(y_proj) <= field_half_size) {
                count_in_field++;
            }
        }
    }

    double fluence = static_cast<double>(count_in_field) / total_particles;
    std::cout << "Fluence at 10x10 cm^2 field (z = 100 cm): " << fluence << std::endl;

    file.close();
    return 0;
}

