echo "BUILD STARTS"
python -m pip install -r requirements.txt
echo "Installing dependencies done."
python3.9 manage.py collectstatic --noinput --clear
echo "Collecting static files done."
echo "BUILD ENDS"
