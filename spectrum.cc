#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
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

int main() {
    std::ifstream file("../Varian_TrueBeam6MV_01.IAEAphsp", std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return 1;
    }

    IAEA_Particle particle;
    char particle_type_raw;
    int total_particles = 0;
    double sum_energy_photon = 0, sum_energy_electron = 0, sum_energy_positron = 0;
    int count_photon = 0, count_electron = 0, count_positron = 0;

    // Histogram data storage
    std::vector<int> photon_hist(100, 0);
    std::vector<int> electron_hist(100, 0);
    std::vector<int> positron_hist(100, 0);
    const float energy_min = 0.0f;
    const float energy_max = 10.0f;
    const float bin_width = 0.1f;

    while (file.read(&particle_type_raw, sizeof(char))) {
        short is = (particle_type_raw < 0) ? -1 : 1;
        particle.particle_type = std::abs(particle_type_raw);

        if (!file.read(reinterpret_cast<char*>(&particle.energy), sizeof(float))) break;
        particle.is_new_history = (particle.energy < 0) ? 1 : 0;
        particle.energy = std::fabs(particle.energy);

        if (!file.read(reinterpret_cast<char*>(&particle.x), sizeof(float) * 5)) break;

        double aux = (particle.u * particle.u + particle.v * particle.v);
        particle.w = (aux <= 1.0) ? static_cast<float>(is * sqrt(1.0 - aux)) : 0.0f;

        if (particle.particle_type == 1) { // Photons
            int bin = static_cast<int>((particle.energy - energy_min) / bin_width);
            if (bin >= 0 && bin < 100) {
                photon_hist[bin]++;
            }
            sum_energy_photon += particle.energy;
            count_photon++;
        } else if (particle.particle_type == 2) { // Electrons
            int bin = static_cast<int>((particle.energy - energy_min) / bin_width);
            if (bin >= 0 && bin < 100) {
                electron_hist[bin]++;
            }
            sum_energy_electron += particle.energy;
            count_electron++;
        } else if (particle.particle_type == 3) { // Positrons
            int bin = static_cast<int>((particle.energy - energy_min) / bin_width);
            if (bin >= 0 && bin < 100) {
                positron_hist[bin]++;
            }
            sum_energy_positron += particle.energy;
            count_positron++;
        }
        total_particles++;
    }

    file.close();

    std::cout << "Total particles read: " << total_particles << "\n";
    std::cout << "\nAverage Energy:\n";
    if (count_photon > 0) std::cout << "Photons: " << sum_energy_photon / count_photon << " MeV\n";
    if (count_electron > 0) std::cout << "Electrons: " << sum_energy_electron / count_electron << " MeV\n";
    if (count_positron > 0) std::cout << "Positrons: " << sum_energy_positron / count_positron << " MeV\n";

    // Save histograms as CSV
    std::ofstream photon_csv("photon_energy_distribution.csv");
    std::ofstream electron_csv("electron_energy_distribution.csv");
    std::ofstream positron_csv("positron_energy_distribution.csv");

    if (photon_csv.is_open()) {
        photon_csv << "Energy (MeV),Frequency\n";
        for (int i = 0; i < 100; ++i) {
            float energy = energy_min + i * bin_width;
            photon_csv << energy << "," << photon_hist[i] << "\n";
        }
        photon_csv.close();
        std::cout << "Photon energy distribution saved to photon_energy_distribution.csv\n";
    }

    if (electron_csv.is_open()) {
        electron_csv << "Energy (MeV),Frequency\n";
        for (int i = 0; i < 100; ++i) {
            float energy = energy_min + i * bin_width;
            electron_csv << energy << "," << electron_hist[i] << "\n";
        }
        electron_csv.close();
        std::cout << "Electron energy distribution saved to electron_energy_distribution.csv\n";
    }

    if (positron_csv.is_open()) {
        positron_csv << "Energy (MeV),Frequency\n";
        for (int i = 0; i < 100; ++i) {
            float energy = energy_min + i * bin_width;
            positron_csv << energy << "," << positron_hist[i] << "\n";
        }
        positron_csv.close();
        std::cout << "Positron energy distribution saved to positron_energy_distribution.csv\n";
    }

    return 0;
}

