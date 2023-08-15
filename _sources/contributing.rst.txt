============================
Bug Reports and Contributing
============================
The BirDePy source code is stored at `github.com/birdepy/birdepy_project <https://github.com/birdepy/birdepy_project>`_ .

To report a bug use the Issues tab on our Github repository. 

.. image:: bugs.png

In order to contribute we recommend following these steps:

1. Fork the birdepy_project repository. 

.. image:: forks.png

2. Clone the forked repo to your local machine by opening the folder where you would like to store the code and running: ::

	git clone https://github.com/<your_username>/birdepy_project

3. Make the desired changes. 

4. Push your changes to your fork with an appropriate comment: ::

	git add .
	git commit -m "fix: write a comment here"
	git push

The pipeline for releasing new versions will only execute if commit 
messages adhere to the correct formatting, as outlined in the `semantic-release documentation <https://semantic-release.gitbook.io/semantic-release/>`_.
For instance, using a commit message like `Fix stuff` or `New features` will not trigger the release pipeline,
whereas properly formatted message such as `fix: description of fix` or `feat: description of feature` necessary for successful execution.
To ensure that changes are accepted for merging into the main branch, they must successfully
pass a set of unit tests. 

5. Open a pull request: 

	.. image:: pulls.png

We strongly prefer that each pull request focuses on a single aspect of the code. 
