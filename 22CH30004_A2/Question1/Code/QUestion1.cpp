#include <bits/stdc++.h>

using namespace std;

class FractionationTower
{
private:
    map<string, double> feed_composition;
    map<string, double> distillate;
    map<string, double> bottoms;
    map<string, double> component_volatility;
    double rel_volatility;

    static double sum_components(const map<string, double> &mixture)
    {
        double total = 0.0;
        for (const auto &entry : mixture)
            total += entry.second;
        return total;
    }

    double calculate_key_volatility(double heavy_key, double light_key) const
    {
        return light_key / heavy_key;
    }

    double solve_underwood_equation(double feed_quality = 1.0) const
    {
        const double feed_total = sum_components(feed_composition);
        vector<double> volatility_ratios, feed_fractions;

        for (const auto &entry : component_volatility)
        {
            volatility_ratios.push_back(entry.second / component_volatility.at("iC5"));
            feed_fractions.push_back(feed_composition.at(entry.first) / feed_total);
        }

        double theta_min = 0.0, theta_max = 2.0;
        const double precision = 1e-6;
        const int max_iterations = 100;

        for (int i = 0; i < max_iterations && (theta_max - theta_min) > precision; ++i)
        {
            double theta_mid = (theta_min + theta_max) / 2.0;
            double equation_sum = 0.0;

            for (size_t idx = 0; idx < volatility_ratios.size(); ++idx)
            {
                equation_sum += (volatility_ratios[idx] * feed_fractions[idx]) / (volatility_ratios[idx] - theta_mid);
            }
            equation_sum -= (1 - feed_quality);

            (equation_sum > 0) ? theta_min = theta_mid : theta_max = theta_mid;
        }
        return (theta_min + theta_max) / 2.0;
    }

public:
    FractionationTower(map<string, double> feed, map<string, double> overhead,
                       map<string, double> residue, map<string, double> volatility)
        : feed_composition(move(feed)), distillate(move(overhead)),
          bottoms(move(residue)), component_volatility(move(volatility))
    {
        rel_volatility = calculate_key_volatility(component_volatility.at("nC4"),
                                                  component_volatility.at("iC5"));
    }

    double compute_fenske_stages() const
    {
        const double distillate_total = sum_components(distillate);
        const double bottoms_total = sum_components(bottoms);

        const double overhead_heavy = distillate.at("iC5") / distillate_total;
        const double overhead_light = distillate.at("nC4") / distillate_total;
        const double residue_heavy = bottoms.at("iC5") / bottoms_total;
        const double residue_light = bottoms.at("nC4") / bottoms_total;

        const double numerator = (overhead_heavy / overhead_light) * (residue_light / residue_heavy);
        return log(numerator) / log(rel_volatility);
    }

    map<string, pair<double, double>> analyze_secondary_components() const
    {
        map<string, pair<double, double>> distribution;
        for (const auto &entry : feed_composition)
        {
            if (entry.first == "iC5" || entry.first == "nC4")
                continue;

            // Check if the component exists in distillate and bottoms
            double distillate_value = (distillate.find(entry.first) != distillate.end()) ? distillate.at(entry.first) : 0.0;
            double bottoms_value = (bottoms.find(entry.first) != bottoms.end()) ? bottoms.at(entry.first) : 0.0;

            distribution[entry.first] = {distillate_value, bottoms_value};
        }
        return distribution;
    }
    double calculate_underwood_reflux() const
    {
        const double theta_param = solve_underwood_equation();
        const double overhead_total = sum_components(distillate);
        double reflux_sum = 0.0;

        for (const auto &entry : distillate)
        {
            const double vol_ratio = component_volatility.at(entry.first) / component_volatility.at("iC5");
            const double frac = entry.second / overhead_total;
            reflux_sum += (vol_ratio * frac) / (vol_ratio - theta_param);
        }
        return reflux_sum - 1;
    }

    double estimate_operating_stages(double min_reflux, double theoretical_stages) const
    {
        const double operating_reflux = 1.3 * min_reflux;
        const double gilliland_param = (operating_reflux - min_reflux) / (operating_reflux + 1);
        const double correlation_factor = 0.28; // From Gilliland chart

        return (theoretical_stages + 1) / (1 - correlation_factor);
    }

    void perform_analysis()
    {
        cout << "Relative volatility between key components: " << rel_volatility << endl;

        const double min_trays = compute_fenske_stages();
        cout << "Minimum theoretical stages required: " << ceil(min_trays) << endl;

        cout << "\nSecondary component distribution analysis:\n";
        auto components = analyze_secondary_components();
        for (const auto &entry : components)
            cout << entry.first << " - Vapor: " << entry.second.first << " | Liquid: " << entry.second.second << endl;

        const double min_reflux = calculate_underwood_reflux();
        cout << "\nMinimum reflux ratio (Underwood): " << min_reflux << endl;

        const double actual_trays = estimate_operating_stages(min_reflux, min_trays);
        cout << "Estimated actual trays required: " << ceil(actual_trays) << endl;
    }
};

int main()
{
    map<string, double> feed = {{"iC4", 12}, {"nC4", 448}, {"iC5", 36}, {"nC5", 15}, {"C6", 23}, {"C7", 39.1}, {"C8", 272.2}, {"C9", 31.0}};
    map<string, double> overhead = {{"iC4", 12}, {"nC4", 442}, {"iC5", 13}, {"nC5", 1}};
    map<string, double> residue = {{"nC4", 6}, {"iC5", 23}, {"nC5", 14}, {"C6", 23}, {"C7", 39.1}, {"C8", 272.2}, {"C9", 31.0}};
    map<string, double> volatility = {{"iC4", 2.5}, {"nC4", 2.1}, {"iC5", 1.0}, {"nC5", 0.83}, {"C6", 0.5}, {"C7", 0.20}, {"C8", 0.10}, {"C9", 0.08}};

    FractionationTower analyzer(feed, overhead, residue, volatility);
    analyzer.perform_analysis();

    return 0;
}