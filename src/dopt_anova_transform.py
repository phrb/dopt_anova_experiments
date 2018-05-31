import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import IntVector, StrVector, BoolVector, Formula, r

base      = importr("base")
utils     = importr("utils")
stats     = importr("stats")
algdesign = importr("AlgDesign")
car       = importr("car")

def isclose(a, b, rel_tol = 1e-09, abs_tol = 0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def opt_federov(design_formula, data, trials):
    output = algdesign.optFederov(Formula(design_formula),
                                  data,
                                  maxIteration = 1000,
                                  nTrials = trials)
    return output

def transform_lm(design, lm_formula):
    print("Power Transform Step:")
    response = lm_formula.split("~")[0].strip()
    variables = lm_formula.split("~")[1].strip()
    r_snippet = """boxcox_t <- powerTransform(%s, data = %s)
    regression <- lm(bcPower(%s, boxcox_t$lambda) ~ %s, data = %s)
    regression""" %(lm_formula, design.r_repr(), response, variables, design.r_repr())
    transformed_lm = robjects.r(r_snippet)

    return transformed_lm

def anova(design, formula):
    regression = stats.lm(Formula(formula), data = design)
    heteroscedasticity_test = car.ncvTest(regression)
    print("Heteroscedasticity Test p-value:")
    print(heteroscedasticity_test.rx("p")[0][0])

    if heteroscedasticity_test.rx("p")[0][0] < 0.05:
        regression = transform_lm(design, formula)
        heteroscedasticity_test = car.ncvTest(regression)
        print("Heteroscedasticity Test p-value:")
        print(heteroscedasticity_test.rx("p")[0][0])

    summary_regression = stats.summary_aov(regression)
    print("Regression Step:")
    print(summary_regression)

    prf_values = {}

    for k, v in zip(base.rownames(summary_regression[0]), summary_regression[0][4]):
        if k.strip() != "Residuals":
            prf_values[k.strip()] = v

    return regression, prf_values

def predict_best(regression, data):
    print("Predicting Best")
    predicted = stats.predict(regression, data)
    predicted_best = predicted.index(min(predicted))

    p_min = min(predicted)
    i = 0
    for k in range(len(predicted)):
        if isclose(predicted[k], p_min, rel_tol = 1e-5):
            i += 1

    print("Identical predictions (tol = 1e-5): {0}".format(i))
    return data.rx(predicted_best, True)

def prune_data(data, predicted_best, fixed_variables):
    print("Pruning Data")
    conditions = []

    for k, v in fixed_variables.items():
        print(predicted_best.rx2(str(k)))
        if conditions == []:
            conditions = data.rx2(str(k)).ro == predicted_best.rx2(str(k))
        else:
            conditions = conditions.ro & (data.rx2(str(k)).ro == predicted_best.rx2(str(k)))

    pruned_data = data.rx(conditions, True)

    print("Dimensions of Pruned Data: " + str(base.dim(pruned_data)).strip())
    return pruned_data

def get_fixed_variables(predicted_best, ordered_prf_keys, fixed_factors, threshold = 2):
    print("Getting fixed variables")
    variables = ordered_prf_keys
    variables = [v.strip("I)(/1 ") for v in variables]

    unique_variables = []

    for v in variables:
        if v not in unique_variables:
            unique_variables.append(v)
        if len(unique_variables) >= threshold:
            break

    fixed_variables = fixed_factors
    for v in unique_variables:
        fixed_variables[v] = predicted_best.rx2(str(v))[0]

    print("Fixed Variables: " + str(fixed_variables))
    return fixed_variables

def prune_model(factors, inverse_factors, ordered_prf_keys, threshold = 2):
    print("Pruning Model")
    variables = ordered_prf_keys
    variables = [v.strip("I)(/1 ") for v in variables]

    unique_variables = []

    for v in variables:
        if v not in unique_variables:
            unique_variables.append(v)
        if len(unique_variables) >= threshold:
            break

    pruned_factors = [f for f in factors if not f in unique_variables]
    pruned_inverse_factors = [f for f in inverse_factors if not f in unique_variables]

    return pruned_factors, pruned_inverse_factors

def dopt_anova_step(response, factors, inverse_factors, data, step_data, fixed_factors, budget):
    full_model     = "".join([" ~ ",
                              " + ".join(factors)])

    if len(inverse_factors) > 0:
        full_model += " + " + " + ".join(["I(1 / {0})".format(f) for f in
            inverse_factors])

    print(full_model)

    design_formula = full_model
    lm_formula     = response[0] + full_model
    trials         = round(2 * (len(factors) + len(inverse_factors) + 1))

    fixed_variables = fixed_factors

    if budget - len(step_data[0]) < 0:
        print("Full data does not fit on budget")
        if trials < len(step_data[0]):
            print("Computing D-Optimal Design")
            output = opt_federov(design_formula, step_data, trials)
            design = output.rx("design")[0]
        else:
            print("Too few data points for a D-Optimal design")
            design = step_data

        used_experiments = len(design[0])
        regression, prf_values = anova(design, lm_formula)
        ordered_prf_keys       = sorted(prf_values, key = prf_values.get)
        predicted_best         = predict_best(regression, step_data)
        fixed_variables        = get_fixed_variables(predicted_best, ordered_prf_keys,
                                                     fixed_factors)
        pruned_data            = prune_data(data, predicted_best, fixed_variables)

        pruned_factors, pruned_inverse_factors = prune_model(factors, inverse_factors,
                                                             ordered_prf_keys)
    else:
        print("Full data fits on budget, picking best value")
        used_experiments = len(step_data[0])
        prf_values = []
        ordered_prf_keys = []
        pruned_data = []
        pruned_factors = []
        pruned_inverse_factors = []
        predicted_best = step_data.rx((step_data.rx2(response[0]).ro == min(step_data.rx(response[0])[0])),
                                  True)

    return {"prf_values": prf_values,
            "ordered_prf_keys": ordered_prf_keys,
            "predicted_best": predicted_best,
            "pruned_data": pruned_data,
            "pruned_factors": pruned_factors,
            "pruned_inverse_factors": pruned_inverse_factors,
            "fixed_factors": fixed_variables,
            "used_experiments": used_experiments}

def dopt_anova():
    data = utils.read_csv("../data/search_space.csv", header = True)

    initial_factors = ["elements_number", "y_component_number",
                       "vector_length", "temporary_size",
                       "load_overlap", "threads_number", "lws_y"]

    initial_inverse_factors = ["y_component_number", "lws_y",
                               "elements_number", "threads_number"]

    response = ["time_per_pixel"]

    data = data.rx(StrVector(initial_factors + response))
    data_best = data.rx((data.rx2(response[0]).ro == min(data.rx(response[0])[0])),
                        True)

    step_factors = initial_factors
    step_inverse_factors = initial_inverse_factors
    step_space = data

    fixed_factors = {}

    initial_budget = 120
    budget = initial_budget
    used_experiments = 0
    iterations = 4

    for i in range(iterations):
        print("Step {0}".format(i))
        if step_space == []:
            break

        step_data = dopt_anova_step(response,
                                    step_factors,
                                    step_inverse_factors,
                                    data,
                                    step_space,
                                    fixed_factors,
                                    budget)

        step_space = step_data["pruned_data"]
        step_factors = step_data["pruned_factors"]
        step_inverse_factors = step_data["pruned_inverse_factors"]
        budget -= step_data["used_experiments"]
        used_experiments += step_data["used_experiments"]
        fixed_factors = step_data["fixed_factors"]

        print("Fixed Factors: " + str(fixed_factors))

        if step_space != []:
            step_best = step_space.rx((step_space.rx2(response[0]).ro ==
                min(step_space.rx(response[0])[0])), True)

            print("Best Step Slowdown: " +
                    str(step_best.rx(response[0])[0][0] /
                        data_best.rx(response[0])[0][0]))

        print("Slowdown: " +
                str(step_data["predicted_best"].rx(response[0])[0][0] /
                    data_best.rx(response[0])[0][0]))
        print("Budget: {0}/{1}".format(used_experiments, initial_budget))

dopt_anova()
