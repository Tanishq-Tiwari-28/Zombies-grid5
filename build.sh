echo "BUILD STARTS"
python -m pip install -r requirements.txt
python manage.py collectstatic --noinput --clear
echo "BUILD ENDS"