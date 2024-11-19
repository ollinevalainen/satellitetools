import nox


@nox.session(python=["3.10", "3.11"])
def test(session):
    session.install(".")
    session.install("pytest")
    session.run("pytest", "tests/")
