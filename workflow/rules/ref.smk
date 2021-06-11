x = 2

def get_E(x):
    E = x + 1
    return E

D = 2

if D == get_E(x):
    rule test_conditional:
        output:
            touch("bla")

if D < get_E(x):
    rule test_conditional2:
        output:
            touch("bla"),
            touch("hidden_touch")


