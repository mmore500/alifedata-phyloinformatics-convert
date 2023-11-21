.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

First Time Contributor?
-----------------------

General information for first-time open source contributors can be found at <https://opensource.guide/how-to-contribute/> and <https://github.com/firstcontributions/first-contributions>.

Filing a bug report or making a feature request — covered elsewhere in the documentation — can also be a great way to contribute to the project.

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/mmore500/alifedata-phyloinformatics-convert/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

alifedata-phyloinformatics-convert could always use more documentation, whether as part of the official alifedata-phyloinformatics-convert docs, in docstrings, or even on the web in blog posts, articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/mmore500/alifedata-phyloinformatics-convert/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `alifedata-phyloinformatics-convert` for local development.

1. Fork the `alifedata-phyloinformatics-convert` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/alifedata-phyloinformatics-convert.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv alifedata-phyloinformatics-convert
    $ cd alifedata-phyloinformatics-convert/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 alifedata-phyloinformatics-convert tests
    $ python setup.py test or pytest
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring, and add the feature to the list in README.rst.
3. Make sure that the continuous integration passes for all supported Python versions at <https://github.com/mmore500/alifedata-phyloinformatics-convert/actions/workflows/CI.yml>.

Tips
----

To run a subset of tests::

$ pytest tests.test_alifedata-phyloinformatics-convert


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed.
Then run::

$ bump2version patch # possible: major / minor / patch
$ git push
$ git push --tags

GitHub Actions will then deploy to PyPI if tests pass.
