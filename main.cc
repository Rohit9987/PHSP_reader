#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

#pragma pack(push, 1)
struct IAEA_Particle {
    short particle_type;
    float energy;
    float x, y, z;
    float u, v, w;
    float statistical_weight;
    char sign_of_w;
    char is_new_history;
};
#pragma pack(pop)

struct Voxel {
    float x, y, z;
    int count = 0;
    double fluence = 0.0;
    double energy_sum = 0.0;
	Voxel(float x_, float y_, float z_) : x(x_), y(y_), z(z_), count(0), fluence(0.0) {};
};

struct Jaw {
    float x_min, x_max;
    float y_min, y_max;
    float z = 100.0f; // 
};

// Simple empirical function for linear attenuation coefficient in water [cm^-1]
float linear_attenuation_water(float energy_MeV) {
    return 0.07f + 0.1f / energy_MeV;
}

int main() {
    std::ifstream file("../Varian_TrueBeam6MV_01.IAEAphsp", std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open PHSP file!" << std::endl;
        return 1;
    }

    std::ifstream voxel_file("voxel_coordinates.txt");
    if (!voxel_file) {
        std::cerr << "Error: Could not open voxel coordinates file!" << std::endl;
        return 1;
    }

    std::vector<Voxel> voxels;
    float vx, vy, vz;
    while (voxel_file >> vx >> vy >> vz) {
        voxels.push_back(Voxel{vx, vy, vz});
    }
    voxel_file.close();

    std::ifstream jaw_file("jaw_settings.txt");
    if (!jaw_file) {
        std::cerr << "Error: Could not open jaw settings file!" << std::endl;
        return 1;
    }

    Jaw jaw;
    jaw_file >> jaw.x_min >> jaw.x_max >> jaw.y_min >> jaw.y_max;
    jaw_file.close();

    const float voxel_size = 0.2f;
    const float half_voxel = voxel_size / 2.0f;

    IAEA_Particle particle;
    char particle_type_raw;
    int total_particles = 0;

    while (file.read(&particle_type_raw, sizeof(char))) {
        short is = (particle_type_raw < 0) ? -1 : 1;
        particle.particle_type = std::abs(particle_type_raw);

        if (!file.read(reinterpret_cast<char*>(&particle.energy), sizeof(float))) break;
        particle.is_new_history = (particle.energy < 0) ? 1 : 0;
        particle.energy = std::fabs(particle.energy);

        if (!file.read(reinterpret_cast<char*>(&particle.x), sizeof(float) * 5)) break;

        double aux = (particle.u * particle.u + particle.v * particle.v);
        particle.w = (aux <= 1.0) ? static_cast<float>(is * sqrt(1.0 - aux)) : 0.0f;

        if (particle.particle_type == 1 && particle.w > 0) {
            float t_jaw = (jaw.z - particle.z) / particle.w;
            float x_at_jaw = particle.x + particle.u * t_jaw;
            float y_at_jaw = particle.y + particle.v * t_jaw;

            // Reject particle if projected to jaw plane is behind source or outside jaw opening
            if (t_jaw < 0) continue;
            if (x_at_jaw < jaw.x_min || x_at_jaw > jaw.x_max ||
                y_at_jaw < jaw.y_min || y_at_jaw > jaw.y_max) {
                continue;
            }

            for (auto& voxel : voxels) {
                float t = (voxel.z - particle.z) / particle.w;
                float x_proj = particle.x + particle.u * t;
                float y_proj = particle.y + particle.v * t;
                float z_proj = voxel.z;

                if (x_proj >= voxel.x - half_voxel && x_proj <= voxel.x + half_voxel &&
                    y_proj >= voxel.y - half_voxel && y_proj <= voxel.y + half_voxel &&
                    z_proj >= voxel.z - half_voxel && z_proj <= voxel.z + half_voxel) {

                    float dx = voxel.x - particle.x;
                    float dy = voxel.y - particle.y;
                    float dz = voxel.z - particle.z;
                    float distance = std::sqrt(dx * dx + dy * dy + dz * dz);

                    float mu = linear_attenuation_water(particle.energy);
                    float weight = std::exp(-mu * distance);

                    voxel.fluence += weight;
                    voxel.energy_sum += particle.energy * weight;
                    voxel.count++;
                }
            }
        }
        total_particles++;
    }
    file.close();

    std::cout << "Total particles read: " << total_particles << "\n";
    for (const auto& voxel : voxels) {
        std::cout << "Voxel at (" << voxel.x << ", " << voxel.y << ", " << voxel.z << ")\n"
                  << "  Photons intersecting: " << voxel.count << "\n"
                  << "  Attenuated fluence: " << voxel.fluence << " [relative units]\n";
        if (voxel.count > 0) {
            std::cout << "  Mean photon energy: " << (voxel.energy_sum / voxel.fluence) << " MeV\n";
        }
        std::cout << "\n";
    }

    return 0;
}

